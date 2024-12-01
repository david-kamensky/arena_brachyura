// Module for low-level pixel operations moving pixels from one surface (e.g., a texture) to a different surface (e.g., the screen).

// TODO: Analogous language features to C++ templates and `constexpr` to avoid branching without copy/paste?

extern crate sdl2;

use sdl2::surface::Surface;
use std::ffi::c_void;

use nalgebra::SVector;

// Most general pixel transfer function; useful for debugging and prototyping, but has too much branching logic
// for optimal use tight loops.  Returns true if the pixel was transferred, false if not (e.g., for transparency
// or out-of-bounds.)
pub fn transfer_pixel_general(source: &Surface, dest: &Surface, sx: i32, sy: i32, dx: i32, dy: i32, transparent: bool) -> bool {

    let source_bpp = source.pixel_format_enum().byte_size_per_pixel() as i32;
    let source_pitch = source.pitch() as i32;
    let source_w = source.width() as i32;
    let source_h = source.height() as i32;

    let dest_bpp = dest.pixel_format_enum().byte_size_per_pixel() as i32;
    let dest_pitch = dest.pitch() as i32;
    let dest_w = dest.width() as i32;
    let dest_h = dest.height() as i32;

    if((dx >= 0) && (dy >=  0) && (dx < dest_w) && (dy < dest_h) && (sx >= 0) && (sy >=  0) && (sx < source_w) && (sy < source_h)){
        let source_offset: i32 = source_bpp*sx + source_pitch*sy;
        let dest_offset: i32 = dest_bpp*dx + dest_pitch*dy;
        unsafe {
            let source_pixels: *const c_void = (*source.raw()).pixels;
            let source_pixels_offset: *const c_void
                = source_pixels.wrapping_add(source_offset as usize);
            if(!(transparent
                 // Use color (0,0,0) to signal transparency.
                 && (*(source_pixels_offset as *const u8)==0)
                 && (*((source_pixels_offset.wrapping_add(1 as usize)) as *const u8)==0)
                 && (*((source_pixels_offset.wrapping_add(2 as usize)) as *const u8)==0))){
                let dest_pixels: *mut c_void = (*dest.raw()).pixels;
                let dest_pixels_offset: *mut c_void = dest_pixels.wrapping_add(dest_offset as usize);
                std::ptr::copy_nonoverlapping(source_pixels_offset, dest_pixels_offset, source_bpp as usize);
                return true;
            } // if
            return false;
        } // unsafe
    } // if
    return false;
}

// Similar to general case, but with no bounds checks.
pub fn transfer_pixel_no_bounds(source: &Surface, dest: &Surface, sx: i32, sy: i32, dx: i32, dy: i32, transparent: bool) -> bool {

    let source_bpp = source.pixel_format_enum().byte_size_per_pixel() as i32;
    let source_pitch = source.pitch() as i32;
    let dest_bpp = dest.pixel_format_enum().byte_size_per_pixel() as i32;
    let dest_pitch = dest.pitch() as i32;
    let source_offset: i32 = source_bpp*sx + source_pitch*sy;
    let dest_offset: i32 = dest_bpp*dx + dest_pitch*dy;
    unsafe {
        let source_pixels: *const c_void = (*source.raw()).pixels;
        let source_pixels_offset: *const c_void
            = source_pixels.wrapping_add(source_offset as usize);
        if(!(transparent
             // Use color (0,0,0) to signal transparency.
             && (*(source_pixels_offset as *const u8)==0)
             && (*((source_pixels_offset.wrapping_add(1 as usize)) as *const u8)==0)
             && (*((source_pixels_offset.wrapping_add(2 as usize)) as *const u8)==0))){
            let dest_pixels: *mut c_void = (*dest.raw()).pixels;
            let dest_pixels_offset: *mut c_void = dest_pixels.wrapping_add(dest_offset as usize);
            std::ptr::copy_nonoverlapping(source_pixels_offset, dest_pixels_offset, source_bpp as usize);
            return true;
        } // if
    } // unsafe
    return false;
}

// Transferring a pixel with no regard for transparency and no bounds checks on surfaces.
pub fn transfer_pixel_opaque(source: &Surface, dest: &Surface, sx: i32, sy: i32, dx: i32, dy: i32) {

    let source_bpp = source.pixel_format_enum().byte_size_per_pixel() as i32;
    let source_pitch = source.pitch() as i32;
    let dest_bpp = dest.pixel_format_enum().byte_size_per_pixel() as i32;
    let dest_pitch = dest.pitch() as i32;
    let source_offset: i32 = source_bpp*sx + source_pitch*sy;
    let dest_offset: i32 = dest_bpp*dx + dest_pitch*dy;
    unsafe {
        let source_pixels: *const c_void = (*source.raw()).pixels;
        let source_pixels_offset: *const c_void = source_pixels.wrapping_add(source_offset as usize);
        let dest_pixels: *mut c_void = (*dest.raw()).pixels;
        let dest_pixels_offset: *mut c_void = dest_pixels.wrapping_add(dest_offset as usize);
        std::ptr::copy_nonoverlapping(source_pixels_offset, dest_pixels_offset, source_bpp as usize);
    } // unsafe
}

// Transferring a pixel with black as transparent, but no bounds checks on surfaces. Returns true if
// pixel transferred, false if not.
pub fn transfer_pixel_transparent(source: &Surface, dest: &Surface, sx: i32, sy: i32, dx: i32, dy: i32) -> bool {

    let source_bpp = source.pixel_format_enum().byte_size_per_pixel() as i32;
    let source_pitch = source.pitch() as i32;
    let dest_bpp = dest.pixel_format_enum().byte_size_per_pixel() as i32;
    let dest_pitch = dest.pitch() as i32;
    let source_offset: i32 = source_bpp*sx + source_pitch*sy;
    let dest_offset: i32 = dest_bpp*dx + dest_pitch*dy;
    unsafe {
        let source_pixels: *const c_void = (*source.raw()).pixels;
        let source_pixels_offset: *const c_void = source_pixels.wrapping_add(source_offset as usize);
        if(!((*(source_pixels_offset as *const u8)==0) &&
             (*((source_pixels_offset.wrapping_add(1 as usize)) as *const u8)==0) &&
             (*((source_pixels_offset.wrapping_add(2 as usize)) as *const u8)==0))){
            let dest_pixels: *mut c_void = (*dest.raw()).pixels;
            let dest_pixels_offset: *mut c_void = dest_pixels.wrapping_add(dest_offset as usize);
            std::ptr::copy_nonoverlapping(source_pixels_offset, dest_pixels_offset, source_bpp as usize);
            return true;
        } // if
    } // unsafe
    return false;
}

// Scaling the brightness of a pixel by a floating-point number; assumed to be in range [0,1] without checking.
pub fn scale_pixel(surf: &Surface, x: i32, y: i32, scale: f32) {

    let bpp = surf.pixel_format_enum().byte_size_per_pixel() as i32;
    let pitch = surf.pitch() as i32;
    let offset: i32 = bpp*x + pitch*y;
    unsafe {
        let pixels: *mut c_void = (*surf.raw()).pixels;
        let pixels_offset: *mut c_void = pixels.wrapping_add(offset as usize);
        let r = *(pixels_offset as *const u8);
        let g = *((pixels_offset.wrapping_add(1 as usize)) as *const u8);
        let b = *((pixels_offset.wrapping_add(2 as usize)) as *const u8);
        let scaled_r = (scale*(r as f32)) as u8;
        let scaled_g = (scale*(g as f32)) as u8;
        let scaled_b = (scale*(b as f32)) as u8;
        std::ptr::copy_nonoverlapping(&scaled_r, pixels_offset as *mut u8, 1);
        std::ptr::copy_nonoverlapping(&scaled_g, pixels_offset.wrapping_add(1 as usize) as *mut u8, 1);
        std::ptr::copy_nonoverlapping(&scaled_b, pixels_offset.wrapping_add(2 as usize) as *mut u8, 1);
    } // unsafe
}

// Brighten each color component additively, given floating-point RGB values. The RGB values are assumed to
// be normalized such that a value of 1.0 corresponds to trying to add 255 to the corresponding pixel color
// component.
pub fn brighten_pixel(surf: &Surface, x: i32, y: i32, rgb: &SVector::<f32,3>) {

    let bpp = surf.pixel_format_enum().byte_size_per_pixel() as i32;
    let pitch = surf.pitch() as i32;
    let offset: i32 = bpp*x + pitch*y;
    unsafe {
        let pixels: *mut c_void = (*surf.raw()).pixels;
        let pixels_offset: *mut c_void = pixels.wrapping_add(offset as usize);
        let mut r_f = *(pixels_offset as *const u8) as f32;
        let mut g_f = *((pixels_offset.wrapping_add(1 as usize)) as *const u8) as f32;
        let mut b_f = *((pixels_offset.wrapping_add(2 as usize)) as *const u8) as f32;
        r_f += 255.0*rgb[0];
        g_f += 255.0*rgb[1];
        b_f += 255.0*rgb[2];
        r_f = f32::min(f32::max(0.0, r_f), 255.0);
        g_f = f32::min(f32::max(0.0, g_f), 255.0);
        b_f = f32::min(f32::max(0.0, b_f), 255.0);
        let r_u = r_f as u8;
        let g_u = g_f as u8;
        let b_u = b_f as u8;
        std::ptr::copy_nonoverlapping(&r_u, pixels_offset as *mut u8, 1);
        std::ptr::copy_nonoverlapping(&g_u, pixels_offset.wrapping_add(1 as usize) as *mut u8, 1);
        std::ptr::copy_nonoverlapping(&b_u, pixels_offset.wrapping_add(2 as usize) as *mut u8, 1);
    } // unsafe
}
