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
use sdl2::rect::Rect;

use nalgebra::SMatrix;
use nalgebra::DMatrix;
use nalgebra::SVector;

use rand::Rng;
use rand::rngs::ThreadRng;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use std::cmp::max;
use std::cmp::min;
use std::ffi::c_void;
use std::vec::Vec;
use std::process::exit;
use std::fs;
use std::env;

// Some global constants:
static W: u32 = 640;
static H: u32 = 480;
// static W: u32 = 800;
// static H: u32 = 600;
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
//static PL_PROJ_SPEED: f32 = 1.7;
static PL_PROJ_SPEED: f32 = 2.0;
static PL_PROJ_MAX_FLIGHT_TIME: i32 = 3000;
static PL_PROJ_R: f32 = 50.0;
static PL_PROJ_SPLASH_TIME: i32 = 250;
static PL_GUN_SCREEN_FRAC: f32 = 0.3;
static PLAYER_ACCEL: f32 = 7e-3;
static PLAYER_JUMP_SPEED: f32 = 0.8;
static GRAV_ACCEL: f32 = 2e-3;
static PLAYER_AIR_CONTROL: f32 = 0.25;
static MOVE_DAMP_TIMESCALE: f32 = 300.0;
static BOB_AMPLITUDE: f32 = 30.0;
static BOB_SPEED: f32 = 0.007;
static MONSTER_SPEED: f32 = 0.6;
static MONSTER_R: f32 = 150.0;
static MONSTER_COLLISION_DAMP: f32 = 0.5;
static MONSTER_COLLISION_ITERS: i32 = 7;
static DAMAGE_IMPULSE: f32 = 2.5;
static DEATH_TIME: i32 = 300;
static SPAWN_TIME: i32 = 300;
static SCORE_CHANGE: f32 = 0.1;
//static SCORE_DECAY_RATE: f32 = 0.00006;
static SCORE_DECAY_RATE: f32 = 0.00004;
static SCORE_START: f32 = 0.5;
static SCORE_BAR_HEIGHT_FRAC: f32 = 0.025;
static SCORE_BAR_MARGIN_FRAC: f32 = 0.05;
static PLAYER_START: SVector<f32,2> = SVector::<f32,2>::new(200.0,200.0);
static DATA_ROOT: &str = "data/";
static DEFAULT_TEXTURE_PATH: &str = "default_textures/";

static PI: f32 = std::f32::consts::PI;

// Derived constants:
static W2: u32 = W/2;
static H2: u32 = H/2;
static RATIO: f32 = (W as f32)/(H as f32);

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
                                                  -(H2 as f32))/(W2 as f32),
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
                         true, false);
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
        if((x0[1] > 0.0 && x0[2] <= RATIO*x0[1]) &&
           (x1[1] > 0.0 && x1[2] <= RATIO*x1[1]) &&
           (x2[1] > 0.0 && x2[2] <= RATIO*x2[1]) &&
           (x3[1] > 0.0 && x3[2] <= RATIO*x3[1])){return;}
        // All vertices below frustum:
        if((x0[1] < 0.0 && x0[2] <= -RATIO*x0[1]) &&
           (x1[1] < 0.0 && x1[2] <= -RATIO*x1[1]) &&
           (x2[1] < 0.0 && x2[2] <= -RATIO*x2[1]) &&
           (x3[1] < 0.0 && x3[2] <= -RATIO*x3[1])){return;}
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

        let j0 = ((W2 as f32)*x0[1]/x0[2]) as i32;
        let j1 = ((W2 as f32)*x1[1]/x1[2]) as i32;
        let j2 = ((W2 as f32)*x2[1]/x2[2]) as i32;
        let j3 = ((W2 as f32)*x3[1]/x3[2]) as i32;
        j_min = max(-(H2 as i32), min(j0,min(j1,min(j2,j3))) - 1);
        j_max = min((H2 as i32), max(j0,max(j1,max(j2,j3))) + 1);
    }

    // Iterate over the screen-space bounding box:

    // Iterate screen rows:
    for j in j_min..j_max {
        A[(1,2)] = ((j as f32)+0.5)/(W2 as f32);
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

// This is for rendering 2D HUD elements and updating the z-buffer to skip
// any 3D geometry behind them.
fn draw_sprite_2d(source: &Surface, dest: &mut Surface, rect: &Rect,
                  z_buffer: &mut DMatrix<f32>){
    let w = rect.width() as f32;
    let h = rect.height() as f32;
    let s_w = source.width() as f32;
    let s_h = source.height() as f32;
    let x = rect.x() as f32;
    let y = rect.y() as f32;
    let i0 = max(x as i32, 0 as i32);
    let i1 = min((x + w) as i32, W as i32);
    let j0 = max(y as i32, 0 as i32);
    let j1 = min((y + h) as i32, H as i32);
    // Iterate over screen pixels covered by `rect` and transfer corresponding
    // sprite pixels.
    for j in j0..j1{
        let s_y = (s_h*((j-(y as i32)) as f32)/h) as i32;
        for i in i0..i1{
            let s_x = (s_w*((i-(x as i32)) as f32)/w) as i32;
            if(transfer_pixel(source, dest, s_x, s_y, i, j, true)){
                // If the pixel is non-transparent, set the z-buffer to zero
                // to block anything from rendering over the 2D sprite.
                z_buffer[(j as usize,i as usize)] = 0.0;
            } // if
        } // i
    } // j
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
    let x_trans = SVector::<f32,3>::new(x[0] - player.x[0], player.z - x[2],
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

enum SpeedLimitType {
    Strict,
    StraferunSoft,
    IsotropicSoft
}

struct Player<'a> {
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
}

impl<'a> Player<'a> {
    fn new(texture_set: &'a TextureSet<'a>) -> Player<'a> {
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
               score_bar_background: &texture_set.score_bar_background,}
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
        let score_bar_rect = Rect::new(0,0,W,
                                       ((H as f32)*3.0*SCORE_BAR_HEIGHT_FRAC)
                                       as u32);
        draw_sprite_2d(self.score_bar_background, dest, &score_bar_rect,
                       z_buffer);
        // Width of gun image on screen:
        let gun_screen_w = PL_GUN_SCREEN_FRAC*(W as f32);
        // Scale height proportionally based on input image:
        let gun_screen_h = (self.gun_sprite.height() as f32)*gun_screen_w
            /(self.gun_sprite.width() as f32);
        let gun_rect = Rect::new((W2 - ((0.5*gun_screen_w) as u32))
                                 .try_into().unwrap(),
                                 (H - (gun_screen_h as u32))
                                 .try_into().unwrap(),
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

struct Monster<'a> {
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

fn collide_monsters_with_player_projectile(monsters: &mut Vec<Monster>,
                                           player: &mut Player){
    if(player.projectile.ready){return;}
    for monster in monsters {
        let monster_x = SVector::<f32,3>::new(monster.x[0], monster.x[1],
                                              monster.z);
        if((monster_x - player.projectile.x).norm() < MONSTER_R){
            player.projectile.reset();
            player.score += SCORE_CHANGE;
            monster.die();
        }
    } // monster
}

struct TextureSet<'a> {
    outer_wall_texture: Surface<'a>,
    inner_wall_texture: Surface<'a>,
    floor_texture: Surface<'a>,
    sky_texture: Surface<'a>,
    monster_sprite: Surface<'a>,
    monster_dead_sprite: Surface<'a>,
    monster_spawn_sprite: Surface<'a>,
    projectile_sprite: Surface<'a>,
    projectile_splash_sprite: Surface<'a>,
    gun_sprite: Surface<'a>,
    gun_ready_sprite: Surface<'a>,
    score_bar_background: Surface<'a>,
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
            score_bar_background: try_override("score_bar_background.png"),}
    }
}

fn parse_texture_list(list_filename: String) -> Vec<String>{
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
    pub inner_wall_texture: &'a Surface<'a>,
    pub outer_wall_texture: &'a Surface<'a>,
    pub floor_texture: &'a Surface<'a>,
    pub sky_texture: &'a Surface<'a>,
}

impl<'a> Level<'a> {
    pub fn new(texture_set: &'a TextureSet<'a>) -> Level<'a> {
        Level{walls: Vec::<Wall>::new(),
              spawns: Vec::<SVector<f32,2>>::new(),
              inner_wall_texture: &texture_set.inner_wall_texture,
              outer_wall_texture: &texture_set.outer_wall_texture,
              floor_texture: &texture_set.floor_texture,
              sky_texture: &texture_set.sky_texture}
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
        let mut last_col: i32 = -1;
        let mut col: i32 = -1;
        let mut factor: f32 = 6.0;
        let WHF = WALL_H as f32;
        // `n_row` x `n_col` grid cells:
        let n_row = 4;
        let n_col = 4;
        let n_pillars_per_row = 2;
        // Add pillars at random columns in each row:
        for row in 0..n_row {
            for pillar in 0..n_pillars_per_row {
                loop {
                    col = rng.gen_range(0..n_col) as i32;
                    if(col != last_col){break;}
                } // loop
                last_col = col;
                self.add_pillar(SVector::<f32,2>::new
                                ((4.0*(col as f32) + 1.0)*WHF,
                                 (4.0*(row as f32) + 1.0)*WHF), rng);
            } // pillar
        } // row

        // Exterior walls:
        for i in 0..(4*n_row) {
            self.add_wall(WallType::ConstantX, (i as f32)*WHF,
                          ((i+1) as f32)*WHF, 0.0, self.outer_wall_texture);
            self.add_wall(WallType::ConstantX, (i as f32)*WHF,
                          ((i+1) as f32)*WHF, 4.0*(n_col as f32)*WHF,
                          self.outer_wall_texture);
        } // i
        for i in 0..(4*n_col) {
            self.add_wall(WallType::ConstantY, (i as f32)*WHF,
                          ((i+1) as f32)*WHF, 0.0, self.outer_wall_texture);
            self.add_wall(WallType::ConstantY, (i as f32)*WHF,
                          ((i+1) as f32)*WHF, 4.0*(n_row as f32)*WHF,
                          self.outer_wall_texture);
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
                          format: PixelFormatEnum)
                          -> Result<Surface<'static>, String> {
    let orig_image = <Surface as LoadSurface>::from_file(filename);
    match orig_image {
        Err(error_message) => return Err(error_message),
        Ok(unwrapped_image) => return unwrapped_image.convert_format(format),
    }
}

enum GameState {
    Quit,
    Continue,
    Win,
    Lose,
}

struct Game<'a> {
    pub monsters: Vec<Monster<'a>>,
    pub player: Player<'a>,
    pub level: Level<'a>,
}

impl<'a> Game<'a> {
    pub fn new(texture_set: &'a TextureSet<'a>,
               rng: &mut SmallRng) -> Game<'a>{
        let mut new_game
            = Game{level: Level::<'a>::new(texture_set),
                   player: Player::<'a>::new(texture_set),
                   monsters: Vec::<Monster<'a>>::new()};
        new_game.level.add_walls_randomly(rng);
        for spawn in &new_game.level.spawns {
            new_game.monsters.push(Monster::new(texture_set, *spawn, rng));
        } // i
        new_game.player.x = PLAYER_START;
        return new_game;
    }

    pub fn update_state(&mut self, dt: i32, event_pump: &mut EventPump,
                        rng: &mut SmallRng)
                        -> GameState {

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
        self.player.score -= SCORE_DECAY_RATE*(dt as f32);
        return GameState::Continue;
    }

    pub fn render(&self, draw_surf: &mut Surface, z_buffer: &mut DMatrix<f32>){
        // Reset z-buffer:
        for i in 0..H{ for j in 0..W{
            z_buffer[(i as usize, j as usize)] = FAR_Z;
        }}
        // Draw all sprites:
        for monster in &self.monsters {
            monster.render(draw_surf, &self.player, z_buffer);
        }
        self.player.render_gun_and_projectile(draw_surf, z_buffer);

        // Render the level's walls:
        self.level.draw_all_walls(draw_surf, z_buffer, &self.player);

        // Draw floor and sky last to skip as much as possible through
        // z-buffer populated by other things.
        transform_and_draw_floor(self.level.floor_texture,
                                 draw_surf, &self.player,
                                 0.0, FLOOR_TILE_W,
                                 0.0, FLOOR_TILE_W, z_buffer);
        transform_and_draw_sky(&self.player, &self.level.sky_texture,
                               draw_surf, z_buffer);
        // Drawing score bar:
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
}

fn data_filename(name: &str) -> String {
    return (DATA_ROOT.to_owned()+name).to_string();
}

struct CmdArgs {
    pub level: u32,
    pub windowed: bool,
}

impl CmdArgs {
    pub fn new() -> CmdArgs {
        let argv: Vec<String> = env::args().collect();
        let argc = argv.len();
        let mut skip_arg = false;

        // Default values of arguments:
        let mut level = 0;
        let mut windowed = false;

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
            }
        } // i
        return CmdArgs{level: level, windowed: windowed};
    }
}

fn main() -> Result<(), String> {

    let cmd_args = CmdArgs::new();
    let start_level = cmd_args.level;
    let fullscreen = !(cmd_args.windowed);

    let mut sdl_context = sdl2::init()?;
    let mut video_subsystem = sdl_context.video()?;

    let mut window = if(fullscreen){
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

    let texture_set_list = parse_texture_list(data_filename("texture_sets"));

    let n_texture_sets = texture_set_list.len();
    let mut texture_set_index = (start_level
                                 % (n_texture_sets as u32)) as usize;
    let mut texture_set = TextureSet::new(
        data_filename(&texture_set_list[texture_set_index]), format);

    let mut level_seed: u64 = start_level as u64;
    let mut rng = SmallRng::seed_from_u64(level_seed);

    let mut game = Game::new(&texture_set, &mut rng);

    let mut z_buffer = DMatrix::<f32>::zeros(H as usize, W as usize);

    let timer = sdl_context.timer()?;
    let mut dt: i32 = 10;
    let mut t: i32 = 0;

    // Main game loop:
    loop {
        let mut t_start: i32 = timer.ticks() as i32;
        let game_state = game.update_state(dt, &mut event_pump, &mut rng);
        match game_state {
            GameState::Quit => break,
            GameState::Win
                => {// Restart with a new random seed after winning:
                    level_seed += 1;
                    drop(game);
                    texture_set_index = ((level_seed as u32)
                                         % (n_texture_sets as u32)) as usize;
                    texture_set = TextureSet::new(
                        data_filename(&texture_set_list[texture_set_index]),
                        format);
                    rng = SmallRng::seed_from_u64(level_seed);
                    game = Game::new(&texture_set, &mut rng);},
            GameState::Lose
                => {// Retry with the same random seed after losing:
                    rng = SmallRng::seed_from_u64(level_seed);
                    game = Game::new(&texture_set, &mut rng);},
            GameState::Continue => {},
        }
        game.render(&mut draw_surf, &mut z_buffer);

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
