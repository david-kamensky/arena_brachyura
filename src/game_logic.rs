extern crate sdl2;

use sdl2::event::Event;
use sdl2::EventPump;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::pixels::PixelFormatEnum;
use sdl2::pixels::PixelFormat;
use sdl2::surface::Surface;
use sdl2::image::LoadSurface;
use sdl2::mouse::MouseButton;
use sdl2::rect::Rect;
use sdl2::render::Canvas;
use sdl2::render::Texture;
use sdl2::video::Window;

use nalgebra::DMatrix;
use nalgebra::SVector;

use rand::Rng;
use rand::rngs::ThreadRng;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use std::cmp::max;
use std::cmp::min;
use std::vec::Vec;
use std::fs;
use std::env;

use crate::constants::*;
use crate::pixel_ops::*;
use crate::rendering::*;

pub struct Projectile<'a> {
    pub sprite: &'a Surface<'a>,
    pub splash_sprite: &'a Surface<'a>,
    pub x: SVector<f32,3>,
    pub v: SVector<f32,3>,
    pub ready: bool,
    pub flight_time: i32,
    pub splash_time: i32,
}

impl<'a> Projectile<'a> {
    pub fn new(sprite: &'a Surface<'a>,
               splash_sprite: &'a Surface<'a>) -> Projectile<'a> {
        Projectile{sprite: sprite, splash_sprite: splash_sprite,
                   x: SVector::<f32,3>::new(0.0,0.0,0.0),
                   v: SVector::<f32,3>::new(0.0,0.0,0.0),
                   ready: true,
                   flight_time: 0, splash_time: 0,}
    }
    pub fn reset(&mut self){
        self.ready = true;
        self.splash_time = PL_PROJ_SPLASH_TIME;
    }
    pub fn advance(self: &mut Projectile<'a>, dt: i32){
        if(self.splash_time < 0){self.splash_time = 0;}
        if(self.splash_time > 0){self.splash_time -= dt;}
        if(self.ready){return;}
        self.flight_time += dt;
        if(self.flight_time > PL_PROJ_MAX_FLIGHT_TIME){self.reset();}
        else{self.x += (dt as f32)*self.v;}
    }
    pub fn collide_with_level(&mut self, level: &Level){
        if(self.ready){return;}
        // Check collision with floor, or flying out into sky.
        if(self.x[2] < (-0.5*(WALL_H as f32) + PL_PROJ_R)){self.reset(); return;}
        if(self.x[2] > (0.5*(WALL_H as f32) - PL_PROJ_R)){
            if(level.has_ceiling){self.reset();}
            return;
        }
        // Check for collisions with walls:
        for wall in &level.walls {
            let self_x_2d = SVector::<f32,2>::new(self.x[0], self.x[1]);
            let x = wall.closest_point(&self_x_2d);
            let orthog = self_x_2d - x;
            let dist = orthog.norm();
            if(dist < PL_PROJ_R){
                self.reset();
                return;
            }
        } // wall
    }
}

pub enum SpeedLimitType {
    Strict,
    StraferunSoft,
    IsotropicSoft
}

pub struct Player<'a> {
    pub x: SVector<f32,2>,
    pub v: SVector<f32,2>,
    // Separated-out for 2.5D physics:
    pub z: f32,
    pub vz: f32,
    pub yaw: f32,
    pub pitch: f32,
    pub l: i32,
    pub r: i32,
    pub u: i32,
    pub d: i32,
    pub want_to_jump: bool,
    pub projectile: Projectile<'a>,
    pub gun_sprite: &'a Surface <'a>,
    pub gun_ready_sprite: &'a Surface <'a>,
    pub score: f32,
    pub score_bar_background: &'a Surface <'a>,
    pub minimap_cover: &'a Surface<'a>,
}

static PITCH_EPS: f32 = 1e-6;
impl<'a> Player<'a> {
    pub fn new(texture_set: &'a TextureSet<'a>) -> Player<'a> {
        Player{x: SVector::<f32,2>::new(0.0,0.0),
               v: SVector::<f32,2>::new(0.0,0.0),
               z: 0.0,
               vz: 0.0,
               yaw: 0.0,
               pitch: 0.0,
               l: 0,
               r: 0,
               u: 0,
               d: 0,
               want_to_jump: false,
               projectile: Projectile::new(&texture_set.projectile_sprite,
                                           &texture_set.
                                           projectile_splash_sprite),
               gun_sprite: &texture_set.gun_sprite,
               gun_ready_sprite: &texture_set.gun_ready_sprite,
               score: SCORE_START,
               score_bar_background: &texture_set.score_bar_background,
               minimap_cover: &texture_set.minimap_cover,}
    }
    pub fn handle_input(self: &mut Player<'a>,
                        event_pump: &mut EventPump) -> bool{
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit{..} |
                Event::KeyDown{keycode: Some(Keycode::Escape), .. }
                => return true,
                Event::MouseMotion{xrel, yrel, .. }
                => {self.yaw += (xrel as f32)*SENSITIVITY;
                    self.pitch = f32::max(-0.5*PI + PITCH_EPS,
                                          f32::min(0.5*PI - PITCH_EPS, self.pitch
                                                   + (yrel as f32)
                                                   *SENSITIVITY))},
                Event::MouseButtonDown{mouse_btn: MouseButton::Left, .. }
                => self.fire_projectile(),
                Event::KeyDown{keycode: Some(Keycode::A), .. } => self.l = 1,
                Event::KeyDown{keycode: Some(Keycode::D), .. } => self.r = 1,
                Event::KeyDown{keycode: Some(Keycode::W), .. } => self.u = 1,
                Event::KeyDown{keycode: Some(Keycode::S), .. } => self.d = 1,
                Event::KeyDown{keycode: Some(Keycode::Space), .. }
                => {self.want_to_jump = true},
                Event::KeyUp{keycode: Some(Keycode::A), .. } => self.l = 0,
                Event::KeyUp{keycode: Some(Keycode::D), .. } => self.r = 0,
                Event::KeyUp{keycode: Some(Keycode::W), .. } => self.u = 0,
                Event::KeyUp{keycode: Some(Keycode::S), .. } => self.d = 0,
                Event::KeyUp{keycode: Some(Keycode::Space), .. }
                => {self.want_to_jump = false},
                _ => {}
            } // match
        } // for
        return false;
    }
    pub fn jump_if_desired(&mut self){
        if(self.want_to_jump && self.z == 0.0){
            self.vz = PLAYER_JUMP_SPEED;
            self.want_to_jump = false;
        }
    }
    pub fn fire_projectile(self: &mut Player<'a>){
        if(!self.projectile.ready){return;}
        self.projectile.x = SVector::<f32,3>::new(self.x[0], self.x[1], self.z);
        let pc = self.pitch.cos();
        let ps = self.pitch.sin();
        let yc = self.yaw.cos();
        let ys = self.yaw.sin();
        self.projectile.v = PL_PROJ_SPEED
            *SVector::<f32,3>::new(pc*ys, pc*yc, -ps);
        self.projectile.x += (NEAR_Z/PL_PROJ_SPEED)*self.projectile.v;
        self.projectile.ready = false;
        self.projectile.flight_time = 0;
        self.projectile.splash_time = 0;
    }
    pub fn render_gun_and_projectile(self: &Player<'a>, dest: &mut Surface,
                                     z_buffer: &mut DMatrix<f32>){
        // Draw the score bar as part of the player HUD.
        let score_bar_rect = Rect::new(0,0,W, ((H as f32)*3.0*SCORE_BAR_HEIGHT_FRAC) as u32);
        draw_sprite_2d(self.score_bar_background, dest, &score_bar_rect, z_buffer);
        // Width of gun image on screen:
        let gun_screen_w = PL_GUN_SCREEN_FRAC*(W as f32);
        // Scale height proportionally based on input image:
        let gun_screen_h = (self.gun_sprite.height() as f32)*gun_screen_w
            /(self.gun_sprite.width() as f32);
        let gun_rect = Rect::new((W2 - ((0.5*gun_screen_w) as u32)).try_into().unwrap(),
                                 (H - (gun_screen_h as u32)).try_into().unwrap(),
                                 gun_screen_w as u32, gun_screen_h as u32);
        if(self.projectile.splash_time > 0){
            transform_and_draw_sprite(self.projectile.splash_sprite, dest,
                                      &self.projectile.x,
                                      2.0*PL_PROJ_R, 2.0*PL_PROJ_R,
                                      self, z_buffer);
        }
        if(self.projectile.ready){
            draw_sprite_2d(self.gun_ready_sprite, dest, &gun_rect, z_buffer);
            return;
        }else{
            draw_sprite_2d(self.gun_sprite, dest, &gun_rect, z_buffer);
        }
        transform_and_draw_sprite(self.projectile.sprite, dest,
                                  &self.projectile.x,
                                  2.0*PL_PROJ_R, 2.0*PL_PROJ_R,
                                  self, z_buffer);
    }
    pub fn enforce_speed_limit(&self, a: &mut SVector::<f32,2>,
                               limit_type: SpeedLimitType){
        match limit_type {
            SpeedLimitType::Strict => {
                // Prevent acceleration from increasing speed:
                let norm_v = self.v.norm();
                if(norm_v > PLAYER_SPEED){
                    let v_hat = self.v/norm_v;
                    let a_dot_v_hat = a.dot(&v_hat);
                    if(a_dot_v_hat > 0.0){
                        *a -= a_dot_v_hat*v_hat;
                    }
                }},
            SpeedLimitType::StraferunSoft => {
                // Allow strafe-running for sqrt(2) speed increase, and
                // building up higher speeds through bunny-hopping:
                let nyaw = -self.yaw;
                let forward = SVector::<f32,2>::new(nyaw.cos(), nyaw.sin());
                let rightward = SVector::<f32,2>::new(-nyaw.sin(), nyaw.cos());
                if(self.v.dot(&forward).abs() > PLAYER_SPEED){
                    *a -= a.dot(&forward)*forward;
                }if(self.v.dot(&rightward).abs() > PLAYER_SPEED){
                    *a -= a.dot(&rightward)*rightward;
                }},
            SpeedLimitType::IsotropicSoft => {
                // Limit asymptotic velocity under constant acceleration
                // in any direction, but still allow speed to increase
                // while changing direction:
                let v_dot_a = self.v.dot(&a);
                if(v_dot_a/a.norm() > PLAYER_SPEED){
                    *a -= v_dot_a*self.v/self.v.norm_squared();
                }},
        } // match
    }
    pub fn update_velocity_and_position(self: &mut Player<'a>, dt: i32){
        let accel_scale = if(self.z == 0.0){1.0}else{PLAYER_AIR_CONTROL};
        let accel = accel_scale*PLAYER_ACCEL;
        let ay = ((self.u-self.d) as f32)*accel;
        let ax = ((self.r-self.l) as f32)*accel;
        let nyaw = -self.yaw;
        let mut a = SVector::<f32,2>::new(ax*nyaw.cos() - ay*nyaw.sin(),
                                          ax*nyaw.sin() + ay*nyaw.cos());

        //self.enforce_speed_limit(&mut a, SpeedLimitType::Strict);
        self.enforce_speed_limit(&mut a, SpeedLimitType::StraferunSoft);
        //self.enforce_speed_limit(&mut a, SpeedLimitType::IsotropicSoft);

        self.jump_if_desired();

        // Apply acceleration:
        self.v += (dt as f32)*a;
        self.vz -= (dt as f32)*GRAV_ACCEL;
        // Implicit exponential integrator for frictional damping:
        if(self.z <= 0.0){
            self.v *= (-(dt as f32)/MOVE_DAMP_TIMESCALE).exp();
        }
        // Integrate position:
        self.x += (dt as f32)*self.v;
        self.z += (dt as f32)*self.vz;

        // Correct for collision with floor:
        if(self.z < 0.0){
            self.z = 0.0;
            self.vz = 0.0;
        }
    }

    pub fn collide_with_wall(self: &mut Player<'a>, wall: &Wall){
        let x = wall.closest_point(&(self.x));
        let orthog = self.x - x;
        let dist = orthog.norm();
        if(dist < PLAYER_R){
            self.x = x + (PLAYER_R/dist)*orthog;
            let orthog_dot_v = orthog.dot(&self.v);
            if(orthog_dot_v < 0.0){
                self.v -= (orthog_dot_v/dist/dist)*orthog;
            }
        }
    }
}

pub struct Monster<'a> {
    pub sprite: &'a Surface<'a>,
    pub dead_sprite: &'a Surface<'a>,
    pub spawn_sprite: &'a Surface<'a>,
    pub x: SVector<f32,2>,
    pub v: SVector<f32,2>,
    pub z: f32, // Separated out from 2D movement physics with `x` and `v`.
    pub x_spawn: SVector<f32,3>,
    pub bob_phase: f32,
    pub dead_timer: i32,
    pub spawn_timer: i32,
}

impl<'a> Monster<'a> {
    pub fn new(texture_set: &'a TextureSet<'a>,
               x: SVector<f32,2>, rng: &mut SmallRng) -> Monster<'a> {
        let bob_phase = rng.gen_range(0..100) as f32;
        Monster{sprite: &texture_set.monster_sprite,
                dead_sprite: &texture_set.monster_dead_sprite,
                spawn_sprite: &texture_set.monster_spawn_sprite,
                x: x, v: SVector::<f32,2>::new(0.0,0.0),
                bob_phase: bob_phase,
                z: BOB_AMPLITUDE*bob_phase.sin(),
                dead_timer: 0,
                spawn_timer: 0,
                x_spawn: SVector::<f32,3>::new(0.0,0.0,0.0)}
    }
    pub fn render(self: &Monster<'a>, dest: &mut Surface,
                  player: &Player, z_buffer: &mut DMatrix<f32>){
        let x_center = SVector::<f32,3>::new(self.x[0], self.x[1], self.z);
        transform_and_draw_sprite(if(self.dead_timer > 0){self.dead_sprite}
                                  else{self.sprite},
                                  dest, &x_center,
                                  2.0*MONSTER_R, 2.0*MONSTER_R,
                                  player, z_buffer);
        if(self.spawn_timer > 0){
            transform_and_draw_sprite(self.spawn_sprite,
                                      dest, &self.x_spawn,
                                      2.0*MONSTER_R, 2.0*MONSTER_R,
                                      player, z_buffer);
        }
    }
    pub fn die(self: &mut Monster<'a>){self.dead_timer = DEATH_TIME;}
    pub fn dead(&self) -> bool {return self.dead_timer > 0;}
    pub fn think_and_move(self: &mut Monster<'a>, target: &mut Player,
                          level: &Level, dt: i32, rng: &mut SmallRng){
        // If dead, countdown until respawning.
        if(self.dead_timer > 0){
            self.dead_timer -= dt;
            if(self.dead_timer <= 0){
                self.dead_timer = 0;
                let spawn_index = rng.gen_range(0..(level.spawns.len()));
                self.x = level.spawns[spawn_index];
                self.spawn_timer = SPAWN_TIME;
                self.x_spawn = SVector::<f32,3>::new(self.x[0],self.x[1],0.0);
            }
            return;
        } // end if dead
        if(self.spawn_timer > 0){
            self.spawn_timer -= dt;
        }
        let x_diff = self.x - target.x;
        let norm_x_diff = x_diff.norm();
        let norm_z_diff = (self.z - target.z).abs();
        let r_sum = PLAYER_R + MONSTER_R;
        if((norm_x_diff < r_sum) && (norm_z_diff < r_sum)){
            self.die();
            let n = x_diff/norm_x_diff;
            target.score -= SCORE_CHANGE;
            target.v -= DAMAGE_IMPULSE*n;
            return;
        } // end if colliding player
        self.bob_phase += BOB_SPEED*(dt as f32);
        self.z = BOB_AMPLITUDE*(1.0 + self.bob_phase.sin());

        // Set trial velocity based on player location:
        for i in 0..2 {
            let dx_i = self.x[i] - target.x[i];
            if(dx_i.abs() > r_sum){
                self.v[i] = -MONSTER_SPEED*(dx_i/dx_i.abs());
            }
        } // i
        self.x += (dt as f32)*self.v;

        for wall in &level.walls {
            let x = wall.closest_point(&(self.x));
            let orthog = self.x - x;
            let dist = orthog.norm();
            if(dist < MONSTER_R){self.x = x + (MONSTER_R/dist)*orthog;}
        } // wall
    }
}

static MONSTER_COLLISION_EPS: f32 = 1e-3;
pub fn collide_monster_pair(monsters: &mut Vec<Monster>, i: usize, j: usize){
    if(monsters[i].dead() || monsters[j].dead()){return;}
    let x_i = monsters[i].x;
    let x_j = monsters[j].x;
    let dx = x_i - x_j;
    let norm_dx = dx.norm();
    if(norm_dx > 2.0*MONSTER_R){return;}
    // If they're too close together, spawn-kill one of them and let it respawn somewhere else;
    // otherwise the small division behaves erratically.
    if(norm_dx < (MONSTER_COLLISION_EPS*MONSTER_R)){monsters[i].die(); return;}
    let x_mid = 0.5*(x_i + x_j);
    let dx_hat = (1.0/norm_dx)*dx;
    let x_i_out = x_mid + MONSTER_R*dx_hat;
    let x_j_out = x_mid - MONSTER_R*dx_hat;
    monsters[i].x += MONSTER_COLLISION_DAMP*(x_i_out - x_i);
    monsters[j].x += MONSTER_COLLISION_DAMP*(x_j_out - x_j);
}

// Brute-force $O(n^2)$ collision detection between monsters:
pub fn collide_monsters_with_each_other(monsters: &mut Vec<Monster>){
    for iteration in 0..MONSTER_COLLISION_ITERS {
        for i in 0..monsters.len(){
            for j in (i+1)..monsters.len(){
                collide_monster_pair(monsters, i, j);
            } // j
        } // i
    } // iteration
}

pub fn collide_monsters_with_player_projectile(monsters: &mut Vec<Monster>,
                                           player: &mut Player){
    if(player.projectile.ready){return;}
    for monster in monsters {
        if(monster.dead()){continue;}
        let monster_x = SVector::<f32,3>::new(monster.x[0], monster.x[1],
                                              monster.z);
        if((monster_x - player.projectile.x).norm() < MONSTER_R){
            player.projectile.reset();
            player.score += SCORE_CHANGE;
            monster.die();
        }
    } // monster
}

pub struct GameParameters {
    pub has_ceiling: bool,
    pub darkening_length_scale: f32,
}

impl GameParameters {
    pub fn new(cfg_filename: String) -> GameParameters {

        // Default values for parameters:
        let mut has_ceiling = false;
        let mut darkening_length_scale = -1.0; // Negative is full-bright.

        println!("Parsing configuration '{}'...", cfg_filename);
        let file_string: String = fs::read_to_string(cfg_filename).unwrap();
        let lines = file_string.split("\n");
        for line in lines {
            if(line.len() == 0){continue;}
            let tokens = line.split(" ");
            let mut i = 0;
            let mut param = "";
            let mut value = "";
            for token in tokens {
                if(token.len() == 0){continue;}
                if(i==0){param = token; i += 1;}
                else if(i==1){value = token; break;}
                else{break;}
            }
            match param {
                "" => {}, // Happens in pure-whitespace lines; ignore.
                "has_ceiling" => {has_ceiling = value.parse::<bool>().unwrap();},
                "darkening_length_scale" => {darkening_length_scale = value.parse::<f32>().unwrap();},
                _ => {println!("  Warning: Unknown parameter '{}'", param);}
            }
        }
        println!("  has_ceiling = {}", has_ceiling);
        println!("  darkening_length_scale = {}", darkening_length_scale);
        return GameParameters{has_ceiling: has_ceiling,
                              darkening_length_scale: darkening_length_scale};
    } // line
}

pub struct TextureSet<'a> {
    pub outer_wall_texture: Surface<'a>,
    pub inner_wall_texture: Surface<'a>,
    pub floor_texture: Surface<'a>,
    pub sky_texture: Surface<'a>,
    pub monster_sprite: Surface<'a>,
    pub monster_dead_sprite: Surface<'a>,
    pub monster_spawn_sprite: Surface<'a>,
    pub projectile_sprite: Surface<'a>,
    pub projectile_splash_sprite: Surface<'a>,
    pub gun_sprite: Surface<'a>,
    pub gun_ready_sprite: Surface<'a>,
    pub score_bar_background: Surface<'a>,
    pub minimap_cover: Surface<'a>
}

impl<'a> TextureSet<'a> {
    pub fn new(directory: String,
               format: PixelFormatEnum) -> TextureSet<'a> {

        // This lambda loads a file from the given directory if it exists,
        // or from a default directory otherwise. Every texture type is
        // assumed to exist in the default set.
        let try_override = |filename: &str| -> Surface<'a> {
            let try_load_texture =  load_image_with_format
                (directory.clone()+"/"+filename, format);
            match try_load_texture {
                Err(e) => {
                    println!("  Using default for '{}'", filename);
                    let default_path = DATA_ROOT.to_owned()
                        + DEFAULT_TEXTURE_PATH + filename;
                    return load_image_with_format(default_path,
                                                  format).unwrap();
                },
                Ok(texture) => {
                    println!("  Overriding default for '{}'", filename);
                    return texture;
                },
            }
        };
        println!("Loading textures from '{}'", directory);
        TextureSet{
            outer_wall_texture: try_override("outer_wall_texture.png"),
            inner_wall_texture: try_override("inner_wall_texture.png"),
            floor_texture: try_override("floor.png"),
            sky_texture: try_override("sky.png"),
            monster_sprite: try_override("monster.png"),
            monster_dead_sprite: try_override("monster_dead_sprite.png"),
            monster_spawn_sprite: try_override("monster_spawn_sprite.png"),
            projectile_sprite: try_override("projectile.png"),
            projectile_splash_sprite: try_override("projectile_splash.png"),
            gun_sprite: try_override("gun.png"),
            gun_ready_sprite: try_override("gun_ready.png"),
            score_bar_background: try_override("score_bar_background.png"),
            minimap_cover: try_override("minimap_cover.png"),}
    }
}

pub fn parse_texture_list(list_filename: String) -> Vec<String>{
    let mut list = Vec::<String>::new();
    let list_file_string: String = fs::read_to_string(list_filename).unwrap();
    let lines = list_file_string.split("\n");
    for line in lines {
        if(line.len() > 0){
            list.push(line.to_string());
            println!("Found texture set '{}'", line);
        }
    }
    return list;
}

pub enum WallType {
    ConstantX,
    ConstantY,
}

struct Wall<'a> {
    pub x1: f32,
    pub x2: f32,
    pub y1: f32,
    pub y2: f32,
    pub wall_type: WallType,
    pub texture: &'a Surface<'a>,
}

impl<'a> Wall<'a> {
    pub fn new(wall_type: WallType, start: f32, end: f32, constant: f32,
               texture: &'a Surface<'a>) -> Wall<'a> {
        match wall_type {
            WallType::ConstantX =>
                Wall{x1: constant, x2: constant, y1: start, y2: end,
                     wall_type: wall_type,
                     texture: texture},
            WallType::ConstantY =>
                Wall{y1: constant, y2: constant, x1: start, x2: end,
                     wall_type: wall_type,
                     texture: texture},
        } // match
    } // new

    pub fn render(self: &Wall<'a>, player: &Player, dest: &mut Surface,
                  z_buffer: &mut DMatrix<f32>){
        transform_and_draw_wall(player, self.x1, self.y1, self.x2, self.y2,
                                self.texture, dest, z_buffer);
    }

    pub fn closest_point(self: &Wall<'a>, x: &SVector<f32,2>)
                         -> SVector<f32,2>{
        let x0 = SVector::<f32,2>::new(self.x1, self.y1);
        let x1 = SVector::<f32,2>::new(self.x2, self.y2);
        let dx = x1 - x0;
        let dx2 = dx.dot(&dx);
        let s = -dx.dot(&(x0 - x))/dx2;
        if(s < 0.0 || s > 1.0){
            let d0 = (x - x0).norm();
            let d1 = (x - x1).norm();
            if(d0 < d1){x0}else{x1}
        }else{
            s*dx + x0
        }
    }
} // impl Wall

struct Level<'a> {
    pub walls: Vec::<Wall<'a>>,
    pub spawns: Vec::<SVector<f32,2>>,
    pub inner_wall_texture: &'a Surface<'a>,
    pub outer_wall_texture: &'a Surface<'a>,
    pub floor_texture: &'a Surface<'a>,
    pub sky_texture: &'a Surface<'a>,
    pub has_ceiling: bool,
}

impl<'a> Level<'a> {
    pub fn new(texture_set: &'a TextureSet<'a>,
               parameters: &'a GameParameters) -> Level<'a> {
        Level{walls: Vec::<Wall>::new(),
              spawns: Vec::<SVector<f32,2>>::new(),
              inner_wall_texture: &texture_set.inner_wall_texture,
              outer_wall_texture: &texture_set.outer_wall_texture,
              floor_texture: &texture_set.floor_texture,
              sky_texture: &texture_set.sky_texture,
              has_ceiling: parameters.has_ceiling}
    }
    pub fn add_wall(self: &mut Level<'a>, wall_type: WallType,
                    start: f32, end: f32,
                    constant: f32, texture: &'a Surface<'a>){
        self.walls.push(Wall::new(wall_type, start, end, constant, texture));
    }

    pub fn add_pillar(self: &mut Level<'a>, x: SVector<f32,2>,
                      rng: &mut SmallRng){
        // Index of random corner to be moved to center:
        let move_index = 2*(rng.gen_range(0..4));
        let mut pts = Vec::<SVector<f32,2>>::new();
        // Generate points in order:
        let WHF = WALL_H as f32;
        pts.push(SVector::<f32,2>::new(x[0]        ,x[1]));
        pts.push(SVector::<f32,2>::new(x[0]        ,x[1]+WHF));
        pts.push(SVector::<f32,2>::new(x[0]        ,x[1]+2.0*WHF));
        pts.push(SVector::<f32,2>::new(x[0]+WHF    ,x[1]+2.0*WHF));
        pts.push(SVector::<f32,2>::new(x[0]+2.0*WHF,x[1]+2.0*WHF));
        pts.push(SVector::<f32,2>::new(x[0]+2.0*WHF,x[1]+WHF));
        pts.push(SVector::<f32,2>::new(x[0]+2.0*WHF,x[1]));
        pts.push(SVector::<f32,2>::new(x[0]+WHF    ,x[1]));
        pts.push(SVector::<f32,2>::new(x[0]        ,x[1]));
        // Move one toward the center and add a spawn point where it used
        // to be:
        let temp_pt = pts[move_index];
        pts[move_index] = SVector::<f32,2>::new(x[0]+WHF, x[1]+WHF);
        let spawn_pt = 0.5*(temp_pt + pts[move_index]);
        self.spawns.push(spawn_pt);
        if((move_index == 0) || (move_index == 8)){
            pts[8] = SVector::<f32,2>::new(x[0]+WHF, x[1]+WHF);
        }
        for i in 1..9 {
            if(pts[i][0] == pts[i-1][0]){
                self.add_wall(WallType::ConstantX, pts[i-1][1], pts[i][1],
                              pts[i][0], self.inner_wall_texture);
            }else{
                self.add_wall(WallType::ConstantY, pts[i-1][0], pts[i][0],
                              pts[i][1], self.inner_wall_texture);
            } // if/else
        } // i
    }

    pub fn add_walls_randomly(self: &mut Level<'a>, rng: &mut SmallRng){
        self.walls = Vec::<Wall>::new();
        self.spawns = Vec::<SVector<f32,2>>::new();
        let mut col: i32 = -1;
        let mut factor: f32 = 6.0;
        let WHF = WALL_H as f32;
        // `n_row` x `n_col` grid cells:
        let n_row = 4;
        let n_col = 4;
        let mut layout = DMatrix::<u8>::zeros(n_row, n_col);
        let n_pillars_per_row = 3; // Must be strictly less than `n_col` to avoid infinite loop.
        // Add pillars at random columns in each row:
        for row in 0..n_row {
            let mut n_pillars = n_pillars_per_row;
            // Avoid potential infinite loop:
            if((row == 0) && (n_pillars == n_col)){n_pillars -= 1;}
            for pillar in 0..n_pillars {
                loop {
                    col = rng.gen_range(0..n_col) as i32;
                    if((layout[(row as usize,col as usize)] == 0)
                       // Prevent adding pillar near origin, right where player spawns, for fair start to game.
                       && (!((row==0) && (col==0)))){break;}
                } // loop
                layout[(row as usize, col as usize)] = 1;
                self.add_pillar(SVector::<f32,2>::new((4.0*(col as f32) + 1.0)*WHF, (4.0*(row as f32) + 1.0)*WHF), rng);
            } // pillar
        } // row

        // Exterior walls:
        for i in 0..(4*n_row) {
            self.add_wall(WallType::ConstantX, (i as f32)*WHF, ((i+1) as f32)*WHF, 0.0, self.outer_wall_texture);
            self.add_wall(WallType::ConstantX, (i as f32)*WHF, ((i+1) as f32)*WHF, 4.0*(n_col as f32)*WHF, self.outer_wall_texture);
        } // i
        for i in 0..(4*n_col) {
            self.add_wall(WallType::ConstantY, (i as f32)*WHF, ((i+1) as f32)*WHF, 0.0, self.outer_wall_texture);
            self.add_wall(WallType::ConstantY, (i as f32)*WHF, ((i+1) as f32)*WHF, 4.0*(n_row as f32)*WHF, self.outer_wall_texture);
        } // i
    }

    pub fn draw_all_walls(&self, dest: &mut Surface, z_buffer: &mut DMatrix<f32>, player: &Player){
        for w in self.walls.iter() {w.render(player, dest, z_buffer);}
    } // draw_all_walls
    pub fn collide_player_with_walls(self: &Level<'a>, player: &mut Player){
        for w in self.walls.iter(){
            player.collide_with_wall(w);
        }
    }
}

pub fn load_image_with_format(filename: String, format: PixelFormatEnum) -> Result<Surface<'static>, String> {
    let orig_image = <Surface as LoadSurface>::from_file(filename);
    match orig_image {
        Err(error_message) => return Err(error_message),
        Ok(unwrapped_image) => return unwrapped_image.convert_format(format),
    }
}

pub enum GameState {
    Starting,
    Quit,
    Continue,
    Win,
    Lose,
}

pub enum MenuChoice {
    Wait,
    Quit,
    Continue,
}

pub struct Menu<'a> {
    pub state: GameState,
    pub image: Surface<'a>,
}

impl<'a> Menu<'a> {
    pub fn new(state: GameState,) -> Menu<'a> {
        let mut filename: String = match state {
            GameState::Starting => data_filename("/menu_screens/start_screen.png"),
            GameState::Win => data_filename("/menu_screens/win_screen.png"),
            GameState::Lose => data_filename("/menu_screens/lose_screen.png"),
            _ => "".to_string(), // NOTE: Should be unreachable; `Menu` only created for above states.
        };
        Menu{state: state, image: <Surface as LoadSurface>::from_file(filename).unwrap()}
    }
    // Press escape to quit, space to continue.
    pub fn process_input(&self, event_pump: &mut EventPump) -> MenuChoice {
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit{..} |
                Event::KeyDown{keycode: Some(Keycode::Escape), .. }
                => return MenuChoice::Quit,
                Event::KeyDown{keycode: Some(Keycode::Space), .. }
                => return MenuChoice::Continue,
                _ => {}
            } // match
        } // for
        return MenuChoice::Wait;
    }
    // Busy-waiting until quit/continue choice is made.
    pub fn wait_loop(&self, event_pump: &mut EventPump, canvas: &mut Canvas<Window>) -> bool {
        let texture_creator = canvas.texture_creator();
        let texture = Texture::from_surface(&self.image, &texture_creator).unwrap();
        canvas.copy(&texture, None, None);
        canvas.present();
        loop {
            match self.process_input(event_pump) {
                MenuChoice::Wait => {},
                MenuChoice::Quit => {return true;},
                MenuChoice::Continue => {return false;}
            }
        }
        // Should be unreachable.
        return true;
    }
}

pub struct Game<'a> {
    pub monsters: Vec<Monster<'a>>,
    pub player: Player<'a>,
    pub level: Level<'a>,
    pub parameters: &'a GameParameters,
    pub cmd_args: &'a CmdArgs,
}

impl<'a> Game<'a> {
    pub fn new(parameters: &'a GameParameters, cmd_args: &'a CmdArgs,
               texture_set: &'a TextureSet<'a>,
               rng: &mut SmallRng) -> Game<'a>{
        let mut new_game
            = Game{parameters: parameters, cmd_args: cmd_args,
                   level: Level::<'a>::new(texture_set, parameters),
                   player: Player::<'a>::new(texture_set),
                   monsters: Vec::<Monster<'a>>::new()};
        new_game.level.add_walls_randomly(rng);
        if(!(cmd_args.no_monsters)){
            for spawn in &new_game.level.spawns {
                new_game.monsters.push(Monster::new(texture_set, *spawn, rng));
            } // i
        }
        new_game.player.x = PLAYER_START;
        return new_game;
    }

    pub fn update_state(&mut self, dt: i32, event_pump: &mut EventPump, rng: &mut SmallRng) -> GameState {

        // Check for quitting, winning, or losing the current game.
        if(self.player.handle_input(event_pump)){return GameState::Quit;}
        if(self.player.score >= 1.0){return GameState::Win;}
        if(self.player.score <= 0.0){return GameState::Lose;}

        // Otherwise, update the player and monster states, and continue.
        self.player.update_velocity_and_position(dt);
        self.level.collide_player_with_walls(&mut self.player);
        for mut monster in &mut self.monsters {
            monster.think_and_move(&mut self.player, &self.level, dt, rng);
        }
        collide_monsters_with_each_other(&mut self.monsters);
        collide_monsters_with_player_projectile(&mut self.monsters,
                                                &mut self.player);
        self.player.projectile.advance(dt);
        self.player.projectile.collide_with_level(&self.level);
        if(!(self.cmd_args.no_score_decay)){
            self.player.score -= SCORE_DECAY_RATE*(dt as f32);
        }
        return GameState::Continue;
    }

    pub fn draw_minimap(&self, draw_surf: &mut Surface, z_buffer: &mut DMatrix<f32>){
        // TODO: Parameterize with constants.
        let minimap_w = 128 as i32;
        let minimap_dot_w = 6;
        let minimap_x = (W as i32) - minimap_w;
        let minimap_y = (H as i32) - minimap_w;
        let minimap_rect = Rect::new(minimap_x, minimap_y, minimap_w as u32, minimap_w as u32);
        draw_surf.fill_rect(minimap_rect, Color::RGB(16,16,16));
        let minimap_center = SVector::<f32,2>::new((minimap_x as f32) + 0.5*(minimap_w as f32),
                                                   (minimap_y as f32) + 0.5*(minimap_w as f32));
        for monster in &self.monsters {
            if(monster.dead()){continue;}
            let x_monst_trans = monster.x - self.player.x;
            let cos_yaw = self.player.yaw.cos();
            let sin_yaw = self.player.yaw.sin();
            let x_monst = SVector::<f32,2>::new(x_monst_trans[0]*cos_yaw - x_monst_trans[1]*sin_yaw,
                                                -(x_monst_trans[0]*sin_yaw + x_monst_trans[1]*cos_yaw));
            let minimap_r = 8.0*(WALL_H as f32);
            let x_monst_screen = 0.5*(minimap_w as f32)*x_monst/minimap_r + minimap_center;
            if(((x_monst_screen[0] as i32) > (minimap_x + minimap_dot_w/2)) &&
               ((x_monst_screen[1] as i32) > (minimap_y + minimap_dot_w/2))){
                let monst_rect = Rect::new((x_monst_screen[0] as i32) - minimap_dot_w/2,
                                           (x_monst_screen[1] as i32) - minimap_dot_w/2,
                                           minimap_dot_w as u32, minimap_dot_w as u32);
                draw_surf.fill_rect(monst_rect, Color::RGB(255,0,0));
            }
        } // monster
        let player_rect = Rect::new((minimap_center[0] as i32) - minimap_dot_w/2,
                                    (minimap_center[1] as i32) - minimap_dot_w/2,
                                    minimap_dot_w as u32, minimap_dot_w as u32);
        draw_surf.fill_rect(player_rect, Color::RGB(0,255,0));

        // Add cover over minimap:
        draw_sprite_2d(&self.player.minimap_cover, draw_surf, &minimap_rect, z_buffer);

        // Fill in z-buffer with zeros, so this can be drawn first to speed up floor rendering; cover
        // image is mostly transparent, so this needs to be done manually.
        for i in minimap_x..(W as i32) {
            for j in minimap_y..(H as i32) {
                z_buffer[(j as usize, i as usize)] = 0.0;
            } // j
        } // i
    }

    pub fn draw_score_bar(&self, draw_surf: &mut Surface){
        let Wf = W as f32;
        let Hf = H as f32;
        let score_rect = Rect::new((Wf*SCORE_BAR_MARGIN_FRAC) as i32,
                                   (Hf*SCORE_BAR_HEIGHT_FRAC) as i32,
                                   (Wf*(1.0-2.0*SCORE_BAR_MARGIN_FRAC)
                                    *self.player.score) as u32,
                                   (Hf*SCORE_BAR_HEIGHT_FRAC) as u32);
        let score_red = (255.0*f32::min(1.0, 2.0*(1.0 - self.player.score)))
            as u8;
        let score_green = (255.0*f32::min(1.0, 2.0*self.player.score)) as u8;
        draw_surf.fill_rect(score_rect, Color::RGB(score_red,score_green,0));
    }

    pub fn render(&self, draw_surf: &mut Surface, z_buffer: &mut DMatrix<f32>){
        // Reset z-buffer:
        for i in 0..H{ for j in 0..W{
            z_buffer[(i as usize, j as usize)] = FAR_Z;
        }}

        // Foreground HUD elements that fill in z-buffer with zeros; draw first to reduce cost:
        self.draw_minimap(draw_surf, z_buffer);
        self.player.render_gun_and_projectile(draw_surf, z_buffer);
        self.draw_score_bar(draw_surf);

        // Draw all sprites:
        for monster in &self.monsters {
            monster.render(draw_surf, &self.player, z_buffer);
        }

        // Render the level's walls:
        self.level.draw_all_walls(draw_surf, z_buffer, &self.player);

        // Draw floor and sky last to skip as much as possible through
        // z-buffer populated by other things.
        transform_and_draw_floor(self.level.floor_texture, draw_surf, &self.player,
                                 0.0, FLOOR_TILE_W, 0.0, FLOOR_TILE_W, z_buffer, false);
        if(self.parameters.has_ceiling){
            transform_and_draw_floor(self.level.sky_texture, draw_surf, &self.player,
                                     0.0, FLOOR_TILE_W, 0.0, FLOOR_TILE_W, z_buffer, true);
        }else{
            transform_and_draw_sky(&self.player, &self.level.sky_texture, draw_surf, z_buffer);
        }

        // Post-processing filter that adds a depth-darkening effect.  Must run last.
        depth_darkening(draw_surf, z_buffer, self.parameters.darkening_length_scale);
    }
}

pub fn data_filename(name: &str) -> String {
    return (DATA_ROOT.to_owned()+name).to_string();
}

pub struct CmdArgs {
    pub level: u32,
    pub windowed: bool,
    pub no_monsters: bool,
    pub no_score_decay: bool,
}

impl CmdArgs {
    pub fn new() -> CmdArgs {
        let argv: Vec<String> = env::args().collect();
        let argc = argv.len();
        let mut skip_arg = false;

        // Default values of arguments:
        let mut level = 0;
        let mut windowed = false;
        let mut no_monsters = false;
        let mut no_score_decay = false;

        for i in 0..argc {
            if(skip_arg){skip_arg = false; continue;}
            if(argv[i] == "--level"){
                if(i<(argc-1)){
                    level = argv[i+1].parse::<u32>()
                        .expect("Invalid number provided for '--level'.");
                    skip_arg = true;
                }
            }else if(argv[i] == "--windowed"){
                windowed = true;
            }else if(argv[i] == "--no_monsters"){
                no_monsters = true;
            }else if(argv[i] == "--no_score_decay"){
                no_score_decay = true;
            }
        } // i
        return CmdArgs{level: level, windowed: windowed, no_monsters: no_monsters, no_score_decay: no_score_decay};
    }
}
