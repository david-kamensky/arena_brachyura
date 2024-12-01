extern crate sdl2;

use sdl2::surface::Surface;
use sdl2::rect::Rect;
use sdl2::video::WindowSurfaceRef;
use sdl2::render::Canvas;

use nalgebra::SMatrix;
use nalgebra::DMatrix;
use nalgebra::SVector;

use std::cmp::max;
use std::cmp::min;

use crate::constants::*;
use crate::game_logic::*;
use crate::pixel_ops::*;

pub struct ScreenState<'a> {
    pub drawing_surface: Surface<'a>,
    pub z_buffer: DMatrix<f32>,
    pub n_buffer: Vec::<SVector::<f32,3>>,
    // FIXME: Use a more memory-efficient representation of the full-bright mask, e.g., indexing into a bit-vec. (Or, use for
    // other per-pixel flags and differentiate w/ bitwise OR?)
    pub bright_mask: DMatrix<u8>,
}

pub fn column_major(i: usize, j: usize) -> usize {
    return j*(H as usize) + i;
}
pub fn row_major(i: usize, j: usize) -> usize {
    return i*(W as usize) + j;
}

impl<'a> ScreenState<'a> {
    pub fn new(screen_surf: &WindowSurfaceRef) -> ScreenState<'a> {
        ScreenState{
            // NOTE: First three bytes of each pixel assumed to be RGB in low-level pixel operations; default for screen surface on only
            // system tested, but may not be guaranteed in general.
            drawing_surface: Surface::new(screen_surf.width(), screen_surf.height(), sdl2::pixels::PixelFormatEnum::RGB24).unwrap(),
            z_buffer: DMatrix::<f32>::zeros(screen_surf.height() as usize, screen_surf.width() as usize),
            n_buffer: vec![SVector::<f32,3>::new(0.0,0.0,0.0); (screen_surf.width()*screen_surf.height()) as usize],
            bright_mask: DMatrix::<u8>::zeros(screen_surf.height() as usize, screen_surf.width() as usize),
        }
    } // new

    pub fn reset(&mut self){
        // Reset z-buffer and full-bright mask:
        for i in 0..H{ for j in 0..W{
            self.z_buffer[(i as usize, j as usize)] = FAR_Z;
            self.bright_mask[(i as usize, j as usize)] = 0;
            self.n_buffer[column_major(i as usize, j as usize)] = SVector::<f32,3>::new(0.0,0.0,0.0);
        }}
    } // reset

    // Read from and write to the normal buffer as if it's logically 2D:
    pub fn n(&self, i: usize, j: usize) -> &SVector::<f32,3>{
        return &(self.n_buffer[column_major(i,j)]);
    } // n
}

// NOTE: This interprets the pixel-space of the sky image as a spherical polar coordinate chart, so objects at higher
// elevation angles become severely distorted.
pub fn transform_and_draw_sky(player: &Player, sky: &Surface, screen_state: &mut ScreenState){

    let z_buffer = &mut screen_state.z_buffer;
    let screen = &mut screen_state.drawing_surface;

    let t_h = sky.height() as f32;
    let t_w = sky.width() as f32;
    let cos_pitch = player.pitch.cos();
    let sin_pitch = player.pitch.sin();
    let cos_yaw = player.yaw.cos();
    let sin_yaw = player.yaw.sin();
    let mut x_screen = SVector::<f32,3>::new(0.0,0.0,1.0);
    let mut x_pitch = SVector::<f32,3>::new(0.0,0.0,0.0);

    let j_max: i32 = min(H as i32, max(0,((H2 as f32) + 0.5 - player.pitch.tan()*(W2 as f32)) as i32));

    //for j in 0..H {
    for j in 0..j_max {
        x_screen[1] = ((j as f32)+0.5 - (H2 as f32))/(W2 as f32);
        x_pitch[1] = x_screen[1]*cos_pitch + x_screen[2]*sin_pitch;

        // Below horizon.
        //if(x_pitch[1] > 0.0){continue;}

        x_pitch[2] = -x_screen[1]*sin_pitch + x_screen[2]*cos_pitch;
        for i in 0..W {
            // Render sky last and skip any pixel that's already covered.
            if(z_buffer[(j as usize, i as usize)] < FAR_Z){continue;}

            x_screen[0] = ((i as f32)+0.5 - (W2 as f32))/(W2 as f32);
            x_pitch[0] = x_screen[0];
            let x_yaw = SVector::<f32,3>::new(x_pitch[0]*cos_yaw + x_pitch[2]*sin_yaw,
                                              x_pitch[1],
                                              -x_pitch[0]*sin_yaw + x_pitch[2]*cos_yaw);
            let x = SVector::<f32,3>::new(x_yaw[0], x_yaw[2], -x_yaw[1]);
            // No sky below horizon.
            //if(x[2] < 0.0){continue;}
            let phi = x[2].atan2((x[0]*x[0] + x[1]*x[1]).sqrt());
            // Positive theta goes to right of player for y-axis pointing out of screen.
            let theta = x[0].atan2(x[1]);
            let tx = t_w*(theta + PI)/(2.0*PI);
            let ty = t_h*(0.5*PI - phi)/(0.5*PI);
            transfer_pixel_opaque(sky, screen, tx as i32, ty as i32, i as i32, j as i32);
        } // i
    } // j
}

pub fn transform_and_draw_panel(player: &Player, panel: &VerticalPanel, screen_state: &mut ScreenState){
    let p0 = SVector::<f32,3>::new(panel.x0[0], panel.x0[1], panel.z_top);
    let p1 = SVector::<f32,3>::new(panel.x1[0], panel.x1[1], panel.z_top);
    let p2 = SVector::<f32,3>::new(panel.x0[0], panel.x0[1], panel.z_bottom);
    let pp0 = transform_for_player(&p0, &player);
    let pp1 = transform_for_player(&p1, &player);
    let pp2 = transform_for_player(&p2, &player);
    render_parallelogram(&pp0, &pp1, &pp2, panel.texture, false, true, false, screen_state);
}

pub fn render_parallelogram(x0: &SVector<f32,3>, x1: &SVector<f32,3>, x2: &SVector<f32,3>, texture: &Surface,
                            full_bright: bool, transparent: bool, tile: bool, screen_state: &mut ScreenState){

    let z_buffer = &mut screen_state.z_buffer;
    let n_buffer = &mut screen_state.n_buffer;
    let bright_mask = &mut screen_state.bright_mask;
    let screen = &mut screen_state.drawing_surface;

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

    // Normal vectors will be normalized after all of them have been set, to avoid
    // unnecessary computations for occluded polygons.
    let n = v1.cross(&v2);

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
            if(transfer_pixel_no_bounds(texture, screen, t_col, t_row, W2_i, H2_j, transparent)){
                let s_row = H2_j as usize;
                let s_col = W2_i as usize;
                z_buffer[(s_row, s_col)] = z;
                bright_mask[(s_row, s_col)] = if(full_bright){1}else{0};
                n_buffer[column_major(s_row, s_col)] = n;
            }
        } // j
    } // i
}

pub fn transform_and_draw_sprite(source: &Surface, x: &SVector<f32,3>, width: f32, height: f32,
                                 player: &Player, full_bright: bool, screen_state: &mut ScreenState){
    let x_trans = transform_for_player(&x, player);
    // Define surface to always face player after transforming center point:
    let x0 = SVector::<f32,3>::new(x_trans[0]-0.5*width,
                                   x_trans[1]-0.5*height, x_trans[2]);
    let x1 = SVector::<f32,3>::new(x_trans[0]+0.5*width,
                                   x_trans[1]-0.5*height, x_trans[2]);
    let x2 = SVector::<f32,3>::new(x_trans[0]-0.5*width,
                                   x_trans[1]+0.5*height, x_trans[2]);
    render_parallelogram(&x0, &x1, &x2, source, full_bright, true, false, screen_state);
}

// This is for rendering 2D HUD elements and updating the z-buffer to skip
// any 3D geometry behind them.
pub fn draw_sprite_2d(source: &Surface, rect: &Rect, screen_state: &mut ScreenState){
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
            if(transfer_pixel_transparent(source, &mut screen_state.drawing_surface, s_x, s_y, i, j)){
                // If the pixel is non-transparent, set the z-buffer to zero
                // to block anything from rendering over the 2D sprite.
                screen_state.z_buffer[(j as usize,i as usize)] = 0.0;
                // All 2D sprites are full-bright.
                screen_state.bright_mask[(j as usize,i as usize)] = 1;
            } // if
        } // i
    } // j
}

pub fn transform_and_draw_floor(source: &Surface, player: &Player, x_low: f32, x_high: f32,
                                y_low: f32, y_high: f32, ceiling: bool, screen_state: &mut ScreenState){
    let z = if(ceiling){0.5*(WALL_H as f32)}else{-0.5*(WALL_H as f32)};
    let x0 = SVector::<f32,3>::new(x_low, y_low, z);
    let x1 = SVector::<f32,3>::new(x_high, y_low, z);
    let x2 = SVector::<f32,3>::new(x_low, y_high, z);
    let p0 = transform_for_player(&x0, &player);
    let p1 = transform_for_player(&x1, &player);
    let p2 = transform_for_player(&x2, &player);
    render_parallelogram(&p0, &p1, &p2, source, false, false, true, screen_state);
}

pub fn transform_for_player(x: &SVector<f32,3>, player: &Player) -> SVector<f32,3> {
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

// Darken `draw_surf` based on `z_buffer`, with brightness drop-off dictated by `length_scale`. Negative length scale
// indicates full-bright lighting, and the function is a no-op.
pub fn depth_darkening(length_scale: f32, screen_state: &mut ScreenState) {
    if(length_scale < 0.0){return;}
    let l2 = length_scale*length_scale;
    for i in 0..H {
        let tan_i = ((i as f32) - (H2 as f32))/(W2 as f32);
        for j in 0..W {
            if(screen_state.bright_mask[(i as usize, j as usize)] != 0){continue;}
            let mut n = screen_state.n_buffer[column_major(i as usize, j as usize)];
            let n2 = n.dot(&n);
            let scale = if(n2 == 0.0){
                0.0 // Just render black if no normal is available.
            }else{
                let tan_j = ((j as f32) - (W2 as f32))/(W2 as f32);
                let z = screen_state.z_buffer[(i as usize, j as usize)];
                let x = SVector::<f32,3>::new(z*tan_j, z*tan_i, z);
                let dist2 = x.dot(&x);

                // Heuristic based on aesthetic considerations; ignoring surface normal looks nicer for
                // low-poly/sprite-based graphics and pixel-art style.
                1.0/((dist2/l2) + 1.0)

                // Cosine scaling of brightness gives very cheap low-poly look here:
                // n /= n2.sqrt();
                // let xn = (x.dot(&n)/dist2.sqrt()).abs();
                // xn/((dist2/l2) + 1.0)
            };
            scale_pixel(&mut screen_state.drawing_surface, j as i32, i as i32, scale);
        } // j
    } // i
}

pub struct DynamicLight {
    pub x: SVector::<f32,3>,
    // *Relative* brightness of RGB components, modulated based on length scale
    pub rgb: SVector::<f32,3>,
    // Intuitive handle to modulate overall strength of light with units of length
    pub length_scale: f32,
}

// TODO: Generalize to iterate over a `&Vec` of dynamic lights, with indices obtained from a `CollisionStructure`.
// TODO: Factor-out some duplicated lines from `depth_darkening`
pub fn apply_dynamic_lighting(light: &DynamicLight, player: &Player, screen_state: &mut ScreenState) {
    let l2 = light.length_scale*light.length_scale;
    let light_x = transform_for_player(&(light.x), player);
    for i in 0..H {
        let tan_i = ((i as f32) - (H2 as f32))/(W2 as f32);
        for j in 0..W {
            if(screen_state.bright_mask[(i as usize, j as usize)] != 0){continue;}
            let mut n = screen_state.n_buffer[column_major(i as usize, j as usize)];
            let n2 = n.dot(&n);
            if(n2 > 0.0){
                let tan_j = ((j as f32) - (W2 as f32))/(W2 as f32);
                let z = screen_state.z_buffer[(i as usize, j as usize)];
                let x = SVector::<f32,3>::new(z*tan_j, z*tan_i, z);
                let dx = x - light_x;
                if(dx.dot(&n)*x.dot(&n) > 0.0){
                    // Only brighten the pixel if the light is on the same side of the polygon as the camera.
                    let dist2 = dx.dot(&dx);

                    // Heuristic based on aesthetic considerations; ignoring surface normal looks nicer for
                    // low-poly/sprite-based graphics and pixel-art style.
                    // let scale = 1.0/((dist2/l2) + 1.0);

                    // Cosine scaling of brightness gives very cheap low-poly look here:
                    n /= n2.sqrt();
                    let xn = (dx.dot(&n)/dist2.sqrt()).abs();
                    let scale = xn/((dist2/l2) + 1.0);

                    brighten_pixel(&mut screen_state.drawing_surface, j as i32, i as i32, &(scale*light.rgb));
                } // if on same side as light
            } // if n available
        } // j
    } // i
} // apply_dynamic_lighting
