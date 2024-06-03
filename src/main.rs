#![allow(warnings)]

extern crate sdl2;

use sdl2::event::Event;
use sdl2::EventPump;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::pixels::PixelFormatEnum;
use sdl2::pixels::PixelFormat;
use sdl2::surface::Surface;
use sdl2::video::WindowSurfaceRef;
use sdl2::video::SwapInterval;
use sdl2::image::LoadSurface;
use sdl2::render::Texture;
use sdl2::mouse::MouseButton;

use nalgebra::SMatrix;
use nalgebra::DMatrix;
use nalgebra::SVector;

use rand::Rng;

use std::cmp::max;
use std::cmp::min;
use std::ffi::c_void;
use std::vec::Vec;
use std::process::exit;

// Some global constants:
static W: u32 = 640;
static H: u32 = 480;
// static W: u32 = 1024;
// static H: u32 = 768;
// static W: u32 = 320;
// static H: u32 = 240;
static FAR_Z: f32 = std::f32::MAX;
static NEAR_Z: f32 = 100.0;
static WALL_H: i32 = 480;
static FLOOR_TILE_W: f32 = (WALL_H as f32);
static PLAYER_SPEED: f32 = 1.2;
static PLAYER_R: f32 = 150.0;
static SENSITIVITY: f32 = 0.003;
static PL_PROJ_SPEED: f32 = 1.7;
static PL_PROJ_MAX_FLIGHT_TIME: i32 = 3000;
static PL_PROJ_R: f32 = 50.0;
static PLAYER_ACCEL: f32 = 7e-3;
static MOVE_DAMP_TIMESCALE: f32 = 300.0;
static BOB_AMPLITUDE: f32 = 30.0;
static BOB_SPEED: f32 = 0.007;
static MONSTER_SPEED: f32 = 0.6;
static MONSTER_R: f32 = 150.0;
static MONSTER_COLLISION_DAMP: f32 = 0.5;
static MONSTER_COLLISION_ITERS: i32 = 7;
static DEATH_TIME: i32 = 300;
static PI: f32 = std::f32::consts::PI;


// Derived constants:
static W2: u32 = W/2;
static H2: u32 = H/2;

// Returns whether or not the pixel was actually transferred.
// (May skip for tranparency or out-of-bounds.)
fn transfer_pixel(source: &Surface, dest: &Surface, sx: i32, sy: i32,
                  dx: i32, dy: i32, transparent: bool) -> bool {

    let source_bpp = source.pixel_format_enum().byte_size_per_pixel() as i32;
    let source_pitch = source.pitch() as i32;
    let source_w = source.width() as i32;
    let source_h = source.height() as i32;

    let dest_bpp = dest.pixel_format_enum().byte_size_per_pixel() as i32;
    let dest_pitch = dest.pitch() as i32;
    let dest_w = dest.width() as i32;
    let dest_h = dest.height() as i32;

    if((dx >= 0) && (dy >=  0) && (dx < dest_w) && (dy < dest_h)
       && (sx >= 0) && (sy >=  0) && (sx < source_w) && (sy < source_h)){
        let source_offset: i32 = source_bpp*sx + source_pitch*sy;
        let dest_offset: i32 = dest_bpp*dx + dest_pitch*dy;
        unsafe {
            let source_pixels: *const c_void = (*source.raw()).pixels;
            let source_pixels_offset: *const c_void
                = source_pixels.wrapping_add(source_offset as usize);
            if(!(transparent
                 // Use color (0,0,0) to signal transparency.
                 && (*(source_pixels_offset as *const u8)==0)
                 && (*((source_pixels_offset.wrapping_add(1 as usize))
                       as *const u8)==0)
                 && (*((source_pixels_offset.wrapping_add(2 as usize))
                       as *const u8)==0))){
                let dest_pixels: *mut c_void = (*dest.raw()).pixels;
                let dest_pixels_offset: *mut c_void
                    = dest_pixels.wrapping_add(dest_offset as usize);
                std::ptr::copy_nonoverlapping(source_pixels_offset,
                                              dest_pixels_offset,
                                              source_bpp as usize);
                return true;
            } // if
            return false;
        } // unsafe
    } // if
    return false;
}

// NOTE: This interprets the pixel-space of the sky image as a spherical polar
// coordinate chart, so objects at higher elevation angles become severely
// distorted.
fn transform_and_draw_sky(player: &Player, sky: &Surface,
                          screen: &mut Surface, z_buffer: &DMatrix<f32>){
    let t_h = sky.height() as f32;
    let t_w = sky.width() as f32;
    let cos_pitch = player.pitch.cos();
    let sin_pitch = player.pitch.sin();
    let cos_yaw = player.yaw.cos();
    let sin_yaw = player.yaw.sin();
    for j in 0..H {
        for i in 0..W {
            // Render last and skip any pixel that's already covered.
            if(z_buffer[(j as usize, i as usize)] < FAR_Z){continue;}

            let screen_r = SVector::<f32,3>::new(((i as f32)+0.5
                                                  -(W2 as f32))/(W2 as f32),
                                                 ((j as f32)+0.5
                                                  -(H2 as f32))/(H2 as f32),
                                                 1.0);
            let x_screen = (1.0/screen_r.norm())*screen_r;
            let x_pitch =  SVector::<f32,3>::new(x_screen[0],
                                                 x_screen[1]*cos_pitch
                                                 + x_screen[2]*sin_pitch,
                                                 -x_screen[1]*sin_pitch
                                                 + x_screen[2]*cos_pitch);
            let x_yaw = SVector::<f32,3>::new(x_pitch[0]*cos_yaw
                                              + x_pitch[2]*sin_yaw,
                                              x_pitch[1],
                                              -x_pitch[0]*sin_yaw
                                              + x_pitch[2]*cos_yaw);
            let x = SVector::<f32,3>::new(x_yaw[0], x_yaw[2], -x_yaw[1]);
            // No sky below horizon.
            if(x[2] < 0.0){continue;}
            let phi = x[2].atan2((x[0]*x[0] + x[1]*x[1]).sqrt());
            // Positive theta goes to right of player for y-axis pointing
            // out of screen.
            let theta = x[0].atan2(x[1]);
            let tx = t_w*(theta + PI)/(2.0*PI);
            let ty = t_h*(0.5*PI - phi)/(0.5*PI);
            transfer_pixel(sky, screen, tx as i32, ty as i32,
                           i as i32, j as i32, false);
        } // i
    } // j
}

fn transform_and_draw_wall(player: &Player, x1: f32, y1: f32, x2: f32, y2: f32,
                           texture: &Surface,
                           screen: &mut Surface, z_buffer: &mut DMatrix<f32>){
    let p0 = SVector::<f32,3>::new(x1, y1, 0.5*(WALL_H as f32));
    let p1 = SVector::<f32,3>::new(x2, y2, 0.5*(WALL_H as f32));
    let p2 = SVector::<f32,3>::new(x1, y1, -0.5*(WALL_H as f32));
    let pp0 = transform_for_player(&p0, &player);
    let pp1 = transform_for_player(&p1, &player);
    let pp2 = transform_for_player(&p2, &player);
    render_parallelogram(&pp0, &pp1, &pp2, texture, screen, z_buffer,
                         false, false);
}

static EPS: f32 = 1e-8;
fn render_parallelogram(x0: &SVector<f32,3>, x1: &SVector<f32,3>,
                        x2: &SVector<f32,3>, texture: &Surface,
                        screen: &Surface, z_buffer: &mut DMatrix<f32>,
                        transparent: bool, tile: bool){
    let v1 = x1 - x0;
    let v2 = x2 - x0;
    let x3 = x0 + v1 + v2;

    // Completely skip rendering surfaces in some easily-detected cases:
    if(!tile){
        // All vertices behind camera:
        if((x0[2] <= 0.0) && (x1[2] <= 0.0) && (x2[2] <= 0.0) && (x3[2] <= 0.0))
        {return;}
        // All vertices to right of frustum:
        if((x0[0] > 0.0 && x0[2] <= x0[0]) &&
           (x1[0] > 0.0 && x1[2] <= x1[0]) &&
           (x2[0] > 0.0 && x2[2] <= x2[0]) &&
           (x3[0] > 0.0 && x3[2] <= x3[0])){return;}
        // All vertices to left of frustum:
        if((x0[0] < 0.0 && x0[2] <= -x0[0]) &&
           (x1[0] < 0.0 && x1[2] <= -x1[0]) &&
           (x2[0] < 0.0 && x2[2] <= -x2[0]) &&
           (x3[0] < 0.0 && x3[2] <= -x3[0])){return;}
        // All vertices above frustum:
        if((x0[1] > 0.0 && x0[2] <= x0[1]) &&
           (x1[1] > 0.0 && x1[2] <= x1[1]) &&
           (x2[1] > 0.0 && x2[2] <= x2[1]) &&
           (x3[1] > 0.0 && x3[2] <= x3[1])){return;}
        // All vertices below frustum:
        if((x0[1] < 0.0 && x0[2] <= -x0[1]) &&
           (x1[1] < 0.0 && x1[2] <= -x1[1]) &&
           (x2[1] < 0.0 && x2[2] <= -x2[1]) &&
           (x3[1] < 0.0 && x3[2] <= -x3[1])){return;}
    }

    let mut A = SMatrix::<f32,3,3>::new(-v1[0], -v2[0], 0.0,
                                        -v1[1], -v2[1], 0.0,
                                        -v1[2], -v2[2], 1.0);
    let mut uvt = SVector::<f32,3>::new(0.0, 0.0, 0.0);
    let t_w = texture.width() as f32;
    let t_h = texture.height() as f32;
    let t_w_i = texture.width() as i32;
    let t_h_i = texture.height() as i32;

    // Compute bounding box in screen-space:

    // Resort to brute-force bounds for large polygons that extend behind
    // the camera.
    let mut i_min: i32 = -(W2 as i32);
    let mut i_max: i32 = (W2 as i32);
    let mut j_min: i32 = -(H2 as i32);
    let mut j_max: i32 = (H2 as i32);

    // FIXME: Find a smarter way to bound screen space when some vertices
    // are behind the camera; these formulas only work if all z-coordinates
    // are positive.
    if((!tile) && (x0[2] > 0.0) && (x1[2] > 0.0)
       && (x2[2] > 0.0) && (x3[2] > 0.0)){
        let i0 = ((W2 as f32)*x0[0]/x0[2]) as i32;
        let i1 = ((W2 as f32)*x1[0]/x1[2]) as i32;
        let i2 = ((W2 as f32)*x2[0]/x2[2]) as i32;
        let i3 = ((W2 as f32)*x3[0]/x3[2]) as i32;
        i_min = max(-(W2 as i32), min(i0,min(i1,min(i2,i3))) - 1);
        i_max = min((W2 as i32), max(i0,max(i1,max(i2,i3))) + 1);

        let j0 = ((H2 as f32)*x0[1]/x0[2]) as i32;
        let j1 = ((H2 as f32)*x1[1]/x1[2]) as i32;
        let j2 = ((H2 as f32)*x2[1]/x2[2]) as i32;
        let j3 = ((H2 as f32)*x3[1]/x3[2]) as i32;
        j_min = max(-(H2 as i32), min(j0,min(j1,min(j2,j3))) - 1);
        j_max = min((H2 as i32), max(j0,max(j1,max(j2,j3))) + 1);
    }

    // Iterate over the screen-space bounding box:

    // Iterate screen rows:
    for j in j_min..j_max {
        A[(1,2)] = ((j as f32)+0.5)/(H2 as f32);
        let H2_j = (H2 as i32) + j;
        // Iterate screen columns (consecutive in memory for fixed row):
        for i in i_min..i_max {
            A[(0,2)] = ((i as f32)+0.5)/(W2 as f32);
            if(A.determinant() == 0.0){continue;}
            uvt = A.try_inverse().unwrap()*x0;
            let u = uvt[0];
            let v = uvt[1];
            let t = uvt[2];
            if((!tile) && ((u <= 0.0) || (v <= 0.0)
                           || (u >= 1.0) || (v >= 1.0))){continue;}
            if(t <= 0.0){continue;}
            let W2_i = (W2 as i32) + i;
            let z = x0[2] + u*v1[2] + v*v2[2];
            if(z > z_buffer[(H2_j as usize,W2_i as usize)]){continue;}
            let mut t_col = (u*t_w) as i32;
            let mut t_row = (v*t_h) as i32;
            if(tile){
                t_col = t_col % t_w_i;
                t_row = t_row % t_h_i;
                if(t_col < 0){t_col += t_w_i;}
                if(t_row < 0){t_row += t_h_i;}
            }
            if(transfer_pixel(texture, screen, t_col, t_row, W2_i, H2_j,
                              transparent)){
                z_buffer[(H2_j as usize, W2_i as usize)] = z;
            }
        } // j
    } // i
}

fn transform_and_draw_sprite(source: &Surface, dest: &mut Surface,
                             x: &SVector<f32,3>, width: f32, height: f32,
                             player: &Player, z_buffer: &mut DMatrix<f32>){
    let x_trans = transform_for_player(&x, player);
    // Define surface to always face player after transforming center point:
    let x0 = SVector::<f32,3>::new(x_trans[0]-0.5*width,
                                   x_trans[1]-0.5*height, x_trans[2]);
    let x1 = SVector::<f32,3>::new(x_trans[0]+0.5*width,
                                   x_trans[1]-0.5*height, x_trans[2]);
    let x2 = SVector::<f32,3>::new(x_trans[0]-0.5*width,
                                   x_trans[1]+0.5*height, x_trans[2]);
    render_parallelogram(&x0, &x1, &x2, source, dest, z_buffer, true, false);
}

fn transform_and_draw_floor(source: &Surface, dest: &mut Surface,
                            player: &Player, x_low: f32, x_high: f32,
                            y_low: f32, y_high: f32,
                            z_buffer: &mut DMatrix<f32>){
    let x0 = SVector::<f32,3>::new(x_low, y_low, -0.5*(WALL_H as f32));
    let x1 = SVector::<f32,3>::new(x_high, y_low, -0.5*(WALL_H as f32));
    let x2 = SVector::<f32,3>::new(x_low, y_high, -0.5*(WALL_H as f32));
    let p0 = transform_for_player(&x0, &player);
    let p1 = transform_for_player(&x1, &player);
    let p2 = transform_for_player(&x2, &player);
    render_parallelogram(&p0, &p1, &p2, source, dest, z_buffer, false, true);
}

fn transform_for_player(x: &SVector<f32,3>, player: &Player) -> SVector<f32,3> {
    let x_trans = SVector::<f32,3>::new(x[0] - player.x[0], -x[2],
                                        x[1] - player.x[1]);
    let cos_yaw = player.yaw.cos();
    let sin_yaw = player.yaw.sin();
    let x_yaw = SVector::<f32,3>::new(x_trans[0]*cos_yaw - x_trans[2]*sin_yaw,
                                      x_trans[1],
                                      x_trans[0]*sin_yaw + x_trans[2]*cos_yaw);
    // return x_yaw;
    let cos_pitch = player.pitch.cos();
    let sin_pitch = player.pitch.sin();
    return SVector::<f32,3>::new(x_yaw[0],
                                 x_yaw[1]*cos_pitch - x_yaw[2]*sin_pitch,
                                 x_yaw[1]*sin_pitch + x_yaw[2]*cos_pitch);
}

struct Projectile<'a> {
    pub sprite: &'a Surface<'a>,
    pub x: SVector<f32,3>,
    pub v: SVector<f32,3>,
    pub ready: bool,
    pub flight_time: i32,
}

impl<'a> Projectile<'a> {
    pub fn new(sprite: &'a Surface<'a>) -> Projectile<'a> {
        Projectile{sprite: sprite,
                   x: SVector::<f32,3>::new(0.0,0.0,0.0),
                   v: SVector::<f32,3>::new(0.0,0.0,0.0),
                   ready: true,
                   flight_time: 0,}
    }
    pub fn reset(&mut self){
        self.ready = true;
    }
    pub fn advance(self: &mut Projectile<'a>, dt: i32){
        if(self.ready){return;}
        self.flight_time += dt;
        if(self.flight_time > PL_PROJ_MAX_FLIGHT_TIME){self.reset();}
        else{self.x += (dt as f32)*self.v;}
    }
    pub fn collide_with_level(&mut self, level: &Level){
        if(self.ready){return;}
        // Check collision with floor, or flying out into sky.
        if(self.x[2] < (-0.5*(WALL_H as f32))){self.reset(); return;}
        if(self.x[2] > (0.5*(WALL_H as f32))){return;}
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

struct Player<'a> {
    pub x: SVector<f32,2>,
    pub v: SVector<f32,2>,
    pub yaw: f32,
    pub pitch: f32,
    pub l: i32,
    pub r: i32,
    pub u: i32,
    pub d: i32,
    pub projectile: Projectile<'a>,
}

impl<'a> Player<'a> {
    fn new(projectile_sprite: &'a Surface<'a>) -> Player<'a> {
        Player{
            x: SVector::<f32,2>::new(0.0,0.0),
            v: SVector::<f32,2>::new(0.0,0.0),
            yaw: 0.0,
            pitch: 0.0,
            l: 0,
            r: 0,
            u: 0,
            d: 0,
            projectile: Projectile::new(projectile_sprite),
        }
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
                    self.pitch = f32::max(-0.5*PI,
                                          f32::min(0.5*PI, self.pitch
                                                   + (yrel as f32)
                                                   *SENSITIVITY))},
                Event::MouseButtonDown{mouse_btn: MouseButton::Left, .. }
                => self.fire_projectile(),
                Event::KeyDown{keycode: Some(Keycode::A), .. } => self.l = 1,
                Event::KeyDown{keycode: Some(Keycode::D), .. } => self.r = 1,
                Event::KeyDown{keycode: Some(Keycode::W), .. } => self.u = 1,
                Event::KeyDown{keycode: Some(Keycode::S), .. } => self.d = 1,
                Event::KeyUp{keycode: Some(Keycode::A), .. } => self.l = 0,
                Event::KeyUp{keycode: Some(Keycode::D), .. } => self.r = 0,
                Event::KeyUp{keycode: Some(Keycode::W), .. } => self.u = 0,
                Event::KeyUp{keycode: Some(Keycode::S), .. } => self.d = 0,
                _ => {}
            } // match
        } // for
        return false;
    }
    pub fn fire_projectile(self: &mut Player<'a>){
        if(!self.projectile.ready){return;}
        self.projectile.x = SVector::<f32,3>::new(self.x[0], self.x[1], 0.0);
        let pc = self.pitch.cos();
        let ps = self.pitch.sin();
        let yc = self.yaw.cos();
        let ys = self.yaw.sin();
        self.projectile.v = PL_PROJ_SPEED
            *SVector::<f32,3>::new(pc*ys, pc*yc, -ps);
        self.projectile.x += (NEAR_Z/PL_PROJ_SPEED)*self.projectile.v;
        self.projectile.ready = false;
        self.projectile.flight_time = 0;
    }
    pub fn render_projectile(self: &mut Player<'a>, dest: &mut Surface,
                             z_buffer: &mut DMatrix<f32>){
        if(self.projectile.ready){return;}
        transform_and_draw_sprite(self.projectile.sprite, dest,
                                  &self.projectile.x,
                                  2.0*PL_PROJ_R, 2.0*PL_PROJ_R,
                                  self, z_buffer);
    }
    pub fn update_velocity_and_position(self: &mut Player<'a>, dt: i32){
        let ay = ((self.u-self.d) as f32)*PLAYER_ACCEL;
        let ax = ((self.r-self.l) as f32)*PLAYER_ACCEL;
        let nyaw = -self.yaw;
        let a = SVector::<f32,2>::new(ax*nyaw.cos() - ay*nyaw.sin(),
                                      ax*nyaw.sin() + ay*nyaw.cos());
        let dir_a = (1.0/a.norm())*a;

        // Unconstrained predictor:
        self.v += (dt as f32)*a;

        // Radial return mapping to enforce speed limit:
        let norm_v = self.v.norm();
        if(norm_v > PLAYER_SPEED){
            self.v *= (PLAYER_SPEED/norm_v);
        }
        // Frictional damping:
        self.v *= (-(dt as f32)/MOVE_DAMP_TIMESCALE).exp();

        // Integrate position:
        self.x += (dt as f32)*self.v;
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

struct Monster<'a> {
    pub sprite: &'a Surface<'a>,
    pub dead_sprite: &'a Surface<'a>,
    pub x: SVector<f32,2>,
    pub v: SVector<f32,2>,
    pub z: f32, // Separated out from 2D movement physics with `x` and `v`.
    pub bob_phase: f32,
    pub dead_timer: i32,
}

impl<'a> Monster<'a> {
    pub fn new(sprite: &'a Surface<'a>, dead_sprite: &'a Surface<'a>,
               x: SVector<f32,2>) -> Monster<'a> {
        let bob_phase = rand::thread_rng().gen_range(0..100) as f32;
        Monster{sprite: sprite, dead_sprite: dead_sprite, x: x,
                v: SVector::<f32,2>::new(0.0,0.0),
                bob_phase: bob_phase,
                z: BOB_AMPLITUDE*bob_phase.sin(),
                dead_timer: 0}
    }
    pub fn render(self: &mut Monster<'a>, dest: &mut Surface,
                  player: &Player, z_buffer: &mut DMatrix<f32>){
        let x_center = SVector::<f32,3>::new(self.x[0], self.x[1], self.z);
        transform_and_draw_sprite(if(self.dead_timer > 0){self.dead_sprite}
                                  else{self.sprite},
                                  dest, &x_center,
                                  2.0*MONSTER_R, 2.0*MONSTER_R,
                                  player, z_buffer);
    }
    pub fn die(self: &mut Monster<'a>){self.dead_timer = DEATH_TIME;}
    pub fn think_and_move(self: &mut Monster<'a>, target: &Player,
                          level: &Level, dt: i32){
        // If dead, countdown until respawning.
        if(self.dead_timer > 0){
            self.dead_timer -= dt;
            if(self.dead_timer <= 0){
                self.dead_timer = 0;
                let spawn_index
                    = rand::thread_rng().gen_range(0..(level.spawns.len()));
                self.x = level.spawns[spawn_index];
            }
            return;
        } // end if dead
        if((self.x - target.x).norm() < PLAYER_R){
            self.die();

            // TODO: Update player health/score.

            return;
        } // end if colliding player
        self.bob_phase += BOB_SPEED*(dt as f32);
        self.z = BOB_AMPLITUDE*(1.0 + self.bob_phase.sin());

        // Set trial velocity based on player location:
        for i in 0..2 {
            let dx_i = self.x[i] - target.x[i];
            if(dx_i.abs() > (PLAYER_R + MONSTER_R)){
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

fn collide_monster_pair(monsters: &mut Vec<Monster>, i: usize, j: usize){
    let x_i = monsters[i].x;
    let x_j = monsters[j].x;
    let dx = x_i - x_j;
    let norm_dx = dx.norm();
    if(norm_dx > 2.0*MONSTER_R){return;}
    let x_mid = 0.5*(x_i + x_j);
    let dx_hat = (1.0/norm_dx)*dx;
    let x_i_out = x_mid + MONSTER_R*dx_hat;
    let x_j_out = x_mid - MONSTER_R*dx_hat;
    monsters[i].x += MONSTER_COLLISION_DAMP*(x_i_out - x_i);
    monsters[j].x += MONSTER_COLLISION_DAMP*(x_j_out - x_j);
}

// Brute-force $O(n^2)$ collision detection between monsters:
fn collide_monsters_with_each_other(monsters: &mut Vec<Monster>){
    for iteration in 0..MONSTER_COLLISION_ITERS {
        for i in 0..monsters.len(){
            for j in (i+1)..monsters.len(){
                collide_monster_pair(monsters, i, j);
            } // j
        } // i
    } // iteration
}

fn collide_monsters_with_projectile(monsters: &mut Vec<Monster>,
                                    projectile: &mut Projectile){
    if(projectile.ready){return;}
    for monster in monsters {
        let monster_x = SVector::<f32,3>::new(monster.x[0], monster.x[1],
                                              monster.z);
        if((monster_x - projectile.x).norm() < MONSTER_R){
            projectile.reset();
            monster.die();
        }
    } // monster
}

enum WallType {
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

    pub fn closest_point(self: &Wall<'a>, x: &SVector<f32,2>) -> SVector<f32,2>{
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
    pub texture1: &'a Surface<'a>,
    pub texture2: &'a Surface<'a>,
    pub floor_texture: &'a Surface<'a>,
    pub sky_texture: &'a Surface<'a>,
}

impl<'a> Level<'a> {
    pub fn new(texture1: &'a Surface<'a>,
               texture2: &'a Surface<'a>,
               floor_texture: &'a Surface<'a>,
               sky_texture: &'a Surface<'a>,) -> Level<'a> {
        Level{walls: Vec::<Wall>::new(),
              spawns: Vec::<SVector<f32,2>>::new(),
              texture1: texture1,
              texture2: texture2,
              floor_texture: floor_texture,
              sky_texture: sky_texture}
    }
    pub fn add_wall(self: &mut Level<'a>, wall_type: WallType,
                    start: f32, end: f32,
                    constant: f32, texture: &'a Surface<'a>){
        self.walls.push(Wall::new(wall_type, start, end, constant, texture));
    }

    pub fn add_pillar(self: &mut Level<'a>, x: SVector<f32,2>){
        // Index of random corner to be moved to center:
        let move_index = 2*(rand::thread_rng().gen_range(0..4));
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
                              pts[i][0], self.texture1);
            }else{
                self.add_wall(WallType::ConstantY, pts[i-1][0], pts[i][0],
                              pts[i][1], self.texture1);
            } // if/else
        } // i
    }

    pub fn add_walls_randomly(self: &mut Level<'a>){
        self.walls = Vec::<Wall>::new();
        self.spawns = Vec::<SVector<f32,2>>::new();
        let mut last_col: i32 = -1;
        let mut col: i32 = -1;
        let mut factor: f32 = 6.0;
        let WHF = WALL_H as f32;
        for i in 0..3 {
            loop {
                col = rand::thread_rng().gen_range(0..4) as i32;
                if(col != last_col){break;}
            } // loop
            last_col = col;
            if(i > 1){factor = 2.0;}
            self.add_pillar(SVector::<f32,2>::new(2.0*(col as f32)*WHF,
                                                  factor*WHF));
        } // i
        for i in 0..9 {
            self.add_wall(WallType::ConstantX, (i as f32)*WHF,
                          ((i+1) as f32)*WHF, 0.0, self.texture2);
            self.add_wall(WallType::ConstantX, (i as f32)*WHF,
                          ((i+1) as f32)*WHF, 9.0*WHF, self.texture2);
            self.add_wall(WallType::ConstantY, (i as f32)*WHF,
                          ((i+1) as f32)*WHF, 0.0, self.texture2);
            self.add_wall(WallType::ConstantY, (i as f32)*WHF,
                          ((i+1) as f32)*WHF, 9.0*WHF, self.texture2);
        } // i
    }

    pub fn draw_all_walls(&self, dest: &mut Surface,
                          z_buffer: &mut DMatrix<f32>,
                          player: &Player){
        for w in self.walls.iter() {w.render(player, dest, z_buffer);}
    } // draw_all_walls
    pub fn collide_player_with_walls(self: &Level<'a>, player: &mut Player){
        for w in self.walls.iter(){
            player.collide_with_wall(w);
        }
    }
}

fn load_image_with_format(filename: String,
                          format: PixelFormatEnum) -> Surface<'static> {
    let orig_image: Surface = <Surface as LoadSurface>::from_file(filename)
        .unwrap();
    return orig_image.convert_format(format).unwrap();
}

static FULLSCREEN: bool = true;
fn main() -> Result<(), String> {
    let mut sdl_context = sdl2::init()?;
    let mut video_subsystem = sdl_context.video()?;

    let mut window = if(FULLSCREEN){
        video_subsystem.window("window", W, H)
            .fullscreen().build().map_err(|e| e.to_string())?
    }else{
        video_subsystem.window("window", W, H)
            .position_centered()
            .input_grabbed()
            .build().map_err(|e| e.to_string())?
    };
    sdl_context.mouse().set_relative_mouse_mode(true);
    sdl_context.mouse().show_cursor(false);

    let mut canvas = window.into_canvas().build().map_err(|e| e.to_string())?;
    let texture_creator = canvas.texture_creator();
    let mut event_pump = sdl_context.event_pump()?;
    let screen_surf = canvas.window_mut().surface(&event_pump)?;
    // `draw_surf` is where everything is rendered. It is then copied to the
    // canvas all at once.  Rendering directly to the Window surface is not
    // robust and leads to screen tearing and/or flickering.
    let mut draw_surf = Surface::new(screen_surf.width(), screen_surf.height(),
                                       screen_surf.pixel_format_enum())?;
    let draw_surf_rect = draw_surf.rect();

    let format = draw_surf.pixel_format_enum();
    let outer_wall_texture: Surface
        = load_image_with_format("data/texture2.bmp".to_string(), format);
    let inner_wall_texture: Surface
        = load_image_with_format("data/texture1.bmp".to_string(), format);
    let floor_texture: Surface
        = load_image_with_format("data/floor.bmp".to_string(), format);
    let sky_texture: Surface
        = load_image_with_format("data/sky.png".to_string(), format);
    let monster_sprite: Surface
        = load_image_with_format("data/monster.bmp".to_string(), format);
    let monster_dead_sprite: Surface
        = load_image_with_format("data/deadimg.bmp".to_string(), format);
    let projectile_sprite: Surface
        = load_image_with_format("data/projimg.bmp".to_string(), format);
    let mut player = Player::new(&projectile_sprite);
    player.x = SVector::<f32,2>::new(200.0, 200.0);

    let mut z_buffer = DMatrix::<f32>::zeros(H as usize, W as usize);

    let timer = sdl_context.timer()?;
    let mut dt: i32 = 10;
    let mut t: i32 = 0;

    let x1: f32 = -0.5*(WALL_H as f32);
    let x2: f32 = 0.5*(WALL_H as f32);
    let y1: f32 = 0.5*(WALL_H as f32);
    let y2: f32 = 1.5*(WALL_H as f32);
    let constant_y = y1;
    let constant_x = x1;

    let mut level = Level::new(&inner_wall_texture, &outer_wall_texture,
                               &floor_texture, &sky_texture);
    level.add_walls_randomly();

    let mut monsters = Vec::<Monster>::new();
    for spawn in &level.spawns {
        monsters.push(Monster::new(&monster_sprite, &monster_dead_sprite,
                                   *spawn));
    } // i

    // Main game loop:
    loop {

        let mut t_start: i32 = timer.ticks() as i32;

        let quitting: bool = player.handle_input(&mut event_pump);
        if(quitting){break;}

        player.update_velocity_and_position(dt);

        level.collide_player_with_walls(&mut player);

        // Reset surface to black:
        draw_surf.fill_rect(draw_surf_rect, Color::RGB(0,0,0));

        // Reset z-buffer:
        for i in 0..H{
            for j in 0..W{
                z_buffer[(i as usize, j as usize)] = FAR_Z;
            }
        }

        // Render the level's walls:
        level.draw_all_walls(&mut draw_surf, &mut z_buffer, &player);

        // Draw the level's floor:
        transform_and_draw_floor(&floor_texture, &mut draw_surf, &player,
                                 0.0, FLOOR_TILE_W,
                                 0.0, FLOOR_TILE_W, &mut z_buffer);

        for monster in &mut monsters {
            monster.think_and_move(&player, &level, dt as i32);
            monster.render(&mut draw_surf, &player, &mut z_buffer);
        } // monster

        collide_monsters_with_each_other(&mut monsters);
        collide_monsters_with_projectile(&mut monsters, &mut player.projectile);

        player.projectile.advance(dt);
        player.projectile.collide_with_level(&level);
        player.render_projectile(&mut draw_surf, &mut z_buffer);

        // Draw sky last:
        transform_and_draw_sky(&player, &sky_texture, &mut draw_surf,
                               &z_buffer);

        // NOTE: This is the only way I've found to get software-rendered
        // images to the screen reliably in both fullscreen and windowed
        // modes, without tearing or flickering artifacts.
        let texture = Texture::from_surface(&draw_surf,
                                            &texture_creator).unwrap();
        canvas.copy(&texture, draw_surf_rect, draw_surf_rect);
        canvas.present();

        dt = (timer.ticks() as i32) - t_start;
        t += dt;

    } // loop

    sdl_context.sdldrop();

    Ok(())
}
