use nalgebra::SVector;

// Some global constants:

// pub static W: u32 = 320;
// pub static H: u32 = 240;
pub static W: u32 = 640;
pub static H: u32 = 480;
// pub static W: u32 = 800;
// pub static H: u32 = 600;
// pub static W: u32 = 1024;
// pub static H: u32 = 768;
pub static FAR_Z: f32 = std::f32::MAX;
pub static NEAR_Z: f32 = 100.0;
pub static WALL_H: i32 = 480;
pub static FLOOR_TILE_W: f32 = (WALL_H as f32);
pub static GRASS_H: f32 = 160.0;
pub static GRASS_PER_WALL: usize = 4;
pub static PLAYER_SPEED: f32 = 1.2;
pub static PLAYER_R: f32 = 150.0;
pub static SENSITIVITY: f32 = 0.003;
pub static PL_PROJ_SPEED: f32 = 2.0;
pub static PL_PROJ_MAX_FLIGHT_TIME: i32 = 3000;
pub static PL_PROJ_R: f32 = 50.0;
pub static PL_PROJ_SPLASH_TIME: i32 = 250;
pub static PL_GUN_SCREEN_FRAC: f32 = 0.3;
pub static PLAYER_ACCEL: f32 = 7e-3;
pub static PLAYER_JUMP_SPEED: f32 = 0.8;
pub static GRAV_ACCEL: f32 = 2e-3;
pub static PLAYER_AIR_CONTROL: f32 = 0.25;
pub static MOVE_DAMP_TIMESCALE: f32 = 300.0;
pub static BOB_AMPLITUDE: f32 = 30.0;
pub static BOB_SPEED: f32 = 0.007;
pub static MONSTER_SPEED: f32 = 0.6;
pub static MONSTER_R: f32 = 150.0;
pub static MONSTER_COLLISION_DAMP: f32 = 0.5;
pub static MONSTER_COLLISION_ITERS: i32 = 7;
pub static DAMAGE_IMPULSE: f32 = 2.5;
pub static DEATH_TIME: i32 = 300;
pub static SPAWN_TIME: i32 = 300;
pub static SCORE_CHANGE: f32 = 0.1;
pub static SCORE_DECAY_RATE: f32 = 0.00004;
pub static SCORE_START: f32 = 0.5;
pub static SCORE_BAR_HEIGHT_FRAC: f32 = 0.025;
pub static SCORE_BAR_MARGIN_FRAC: f32 = 0.05;
pub static PLAYER_START: SVector<f32,2> = SVector::<f32,2>::new(200.0,200.0);
pub static DATA_ROOT: &str = "data/";
pub static DEFAULT_TEXTURE_PATH: &str = "default_textures/";
pub static PARAM_FILENAME: &str = "params.cfg";
pub static INITIAL_DT: i32 = 10;
pub static Z_EPS: f32 = 1e-6;

pub static PI: f32 = std::f32::consts::PI;

// Derived constants:
pub static W2: u32 = W/2;
pub static H2: u32 = H/2;
pub static RATIO: f32 = (W as f32)/(H as f32);
