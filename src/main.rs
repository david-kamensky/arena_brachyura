#![allow(warnings)]

extern crate sdl2;

use sdl2::surface::Surface;
use sdl2::video::WindowSurfaceRef;
use sdl2::render::Texture;
use sdl2::mouse::MouseButton;
use sdl2::rect::Rect;

use nalgebra::DMatrix;

use rand::Rng;
use rand::rngs::ThreadRng;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use std::cmp::max;
use std::cmp::min;
use std::vec::Vec;
use std::process::exit;

// Modules:
mod constants;
mod pixel_ops;
mod rendering;
mod game_logic;

use crate::constants::*;
use crate::game_logic::*;

fn main() -> Result<(), String> {

    let cmd_args = CmdArgs::new();
    let start_level = cmd_args.level;
    let fullscreen = !(cmd_args.windowed);

    let mut sdl_context = sdl2::init()?;
    let mut video_subsystem = sdl_context.video()?;

    let mut window = if(fullscreen){
        video_subsystem.window("window", W, H).fullscreen().build().map_err(|e| e.to_string())?
    }else{
        video_subsystem.window("window", W, H).position_centered().input_grabbed().build().map_err(|e| e.to_string())?
    };
    sdl_context.mouse().set_relative_mouse_mode(true);
    sdl_context.mouse().show_cursor(false);

    let mut canvas = window.into_canvas().build().map_err(|e| e.to_string())?;
    let texture_creator = canvas.texture_creator();
    let mut event_pump = sdl_context.event_pump()?;
    let screen_surf = canvas.window_mut().surface(&event_pump)?;
    // `draw_surf` is where everything is rendered. It is then copied to the canvas all at once.  Rendering directly to the Window surface is not
    // robust and leads to screen tearing and/or flickering.
    let mut draw_surf = Surface::new(screen_surf.width(), screen_surf.height(), screen_surf.pixel_format_enum())?;
    let draw_surf_rect = draw_surf.rect();

    let format = draw_surf.pixel_format_enum();

    let texture_set_list = parse_texture_list(data_filename("texture_sets"));

    let n_texture_sets = texture_set_list.len();
    let mut texture_set_index = (start_level % (n_texture_sets as u32)) as usize;
    let mut level_directory = &texture_set_list[texture_set_index];
    let mut texture_set = TextureSet::new(data_filename(level_directory), format);
    let mut cfg_filename = level_directory.to_owned() + "/" + PARAM_FILENAME;
    let mut parameters = GameParameters::new(data_filename(&cfg_filename));

    let mut level_seed: u64 = start_level as u64;
    let mut rng = SmallRng::seed_from_u64(level_seed);

    let mut game = Game::new(&parameters, &cmd_args, &texture_set, &mut rng);

    let mut z_buffer = DMatrix::<f32>::zeros(H as usize, W as usize);

    let timer = sdl_context.timer()?;
    let mut dt: i32 = INITIAL_DT;

    let menu = Menu::new(GameState::Starting);
    if(menu.wait_loop(&mut event_pump, &mut canvas)){
        sdl_context.sdldrop();
        return Ok(());
    }

    // Main game loop:
    loop {
        let mut t_start: i32 = timer.ticks() as i32;
        let game_state = game.update_state(dt, &mut event_pump, &mut rng);
        match game_state {
            GameState::Starting => {},
            GameState::Quit => break,
            GameState::Win
                => {let menu = Menu::new(game_state);
                    if(menu.wait_loop(&mut event_pump, &mut canvas)){break;}
                    // Restart with a new random seed after winning:
                    level_seed += 1;
                    t_start = timer.ticks() as i32;
                    dt = INITIAL_DT;
                    drop(game);
                    texture_set_index = ((level_seed as u32) % (n_texture_sets as u32)) as usize;
                    level_directory = &texture_set_list[texture_set_index];
                    texture_set = TextureSet::new(data_filename(level_directory), format);
                    cfg_filename = level_directory.to_owned() + "/" + PARAM_FILENAME;
                    parameters = GameParameters::new(data_filename(&cfg_filename));
                    rng = SmallRng::seed_from_u64(level_seed);
                    game = Game::new(&parameters, &cmd_args, &texture_set, &mut rng);},
            GameState::Lose
                => {let menu = Menu::new(game_state);
                    if(menu.wait_loop(&mut event_pump, &mut canvas)){break;}
                    // Retry with the same random seed after losing:
                    t_start = timer.ticks() as i32;
                    dt = INITIAL_DT;
                    rng = SmallRng::seed_from_u64(level_seed);
                    game = Game::new(&parameters, &cmd_args, &texture_set, &mut rng);},
            GameState::Continue => {},
        }
        game.render(&mut draw_surf, &mut z_buffer);

        // NOTE: This is the only way I've found to get software-rendered images to the screen reliably in both fullscreen and windowed
        // modes, without tearing or flickering artifacts.
        let texture = Texture::from_surface(&draw_surf, &texture_creator).unwrap();
        canvas.copy(&texture, draw_surf_rect, draw_surf_rect);
        canvas.present();

        dt = (timer.ticks() as i32) - t_start;
    } // loop

    sdl_context.sdldrop();

    Ok(())
}
