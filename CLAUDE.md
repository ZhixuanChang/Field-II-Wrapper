# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

There are general laws should be strictly obeyed:
1. Think before acting. Read existing files before writing code.
2. Be concise in output but thorough in reasoning.
3. Prefer editing over rewriting whole files.
4. Do not re-read files you have already read.
5. Test your code before declaring done.
6. No sycophantic openers or closing fluff.
7. Keep solutions simple and direct.
8. User instructions always override this file.

## Project Overview

Field-II-Wrapper is a MATLAB wrapper library for Field-II (https://field-ii.dk/), designed to simplify 3D ultrasound simulation in uniform media. Field-II uses spatial impulse response methods to simulate ultrasound scattering and acoustic field propagation with point scatterers under Born approximation (no multiple scattering).

The wrapper provides high-level functions for creating complex transducer geometries and managing simulation workflows.

## Core Architecture

### Aperture Generation Functions (xdc_*)

The library follows a hierarchical design for transducer aperture creation:

**Low-level building blocks:**
- `xdc_single_rect.m` - Creates a single rectangular physical element subdivided into mathematical elements
- `xdc_concatenate.m` - Combines multiple apertures into one by merging physical elements

**High-level array generators:**
- `xdc_matrix_array.m` - 2D matrix arrays aligned with x/y axes
- `xdc_sparse_array_rect.m` - Sparse arrays with arbitrary element positions
- `xdc_rect_array.m` - General rectangular arrays with arbitrary element positions
- `xdc_rca_array.m` - Row-Column Addressed (RCA) arrays (incomplete)

All aperture functions call Field-II's native `xdc_rectangles()` function, which requires:
1. `rect` matrix - mathematical element specifications (corners, dimensions, physical element ID)
2. `center` matrix - physical element centers
3. `focus` - fixed focal point [x, y, z]

**Key architectural constraint:** Current implementation requires rectangular elements aligned with x/y axes. The readme.md mentions plans to support rotated arrays using `xdc_triangles` interface with automatic mesh generation.

### Signal Processing Functions

**RF data combination:**
- `f_rf_comb.m` - Aligns RF data from multiple emission cycles by compensating for varying start times and zero-padding to create uniform 3D matrix [samples × elements × emissions]

**Full simulation wrapper:**
- `calc_rf.m` - High-level function that handles complete simulation workflow:
  - Configures Field-II parameters (fs, excitation, impulse response)
  - Loops through transmit events with different delays/apodization
  - Calls `calc_scat_multi()` for each emission
  - Aligns all RF data using time synchronization

### Coordinate System

- X-axis: Element width direction (for linear/matrix arrays)
- Y-axis: Element height direction
- Z-axis: Acoustic propagation direction (normal to transducer face)
- Arrays are centered at origin by default unless `center` parameter is specified

### Mathematical vs Physical Elements

Field-II requires physical elements to be subdivided into smaller "mathematical elements" for accurate simulation:
- Physical elements: Actual transducer elements (e.g., 64 elements in linear array)
- Mathematical elements: Computational subdivisions (e.g., 5×10 sub-elements per physical element)
- The `n_sub_x` and `n_sub_y` parameters control subdivision density
- Aperture functions handle this subdivision internally

## Development Workflow

### Running Examples

Examples demonstrate different array types and imaging modes:
- `example_linear_array_pwi.m` - Plane wave imaging with linear array
- `example_sparse_array_rect.m` - Sparse array 3D imaging
- `example_rca_array_sector_scan.m` - RCA array sector scanning
- `example_concatenate*.m` - Aperture concatenation demonstrations

Run in MATLAB:
```matlab
field_init  % Initialize Field-II (must be in path)
example_name
```

### Typical Simulation Workflow

1. **Initialize Field-II:**
   ```matlab
   field_init
   set_field('fs', fs);  % Sampling frequency
   set_field('c', c);    % Sound speed
   ```

2. **Create apertures:**
   ```matlab
   tx_aper = xdc_matrix_array(nx, pitch_x, width_x, ny, pitch_y, width_y, n_sub_x, n_sub_y, focus);
   rx_aper = xdc_matrix_array(...);
   ```

3. **Configure signals:**
   ```matlab
   xdc_impulse(tx_aper, impulse_response);
   xdc_excitation(tx_aper, excitation);
   ```

4. **Run simulation:**
   ```matlab
   % Option A: Use calc_rf wrapper
   [rfdata, t_start] = calc_rf(fc, fs, tx_aper, rx_aper, tx_delay, tx_apod, scat_pos, scat_amp);

   % Option B: Manual loop with f_rf_comb
   for i = 1:num_emissions
       xdc_apodization(tx_aper, 0, tx_apod(i,:));
       xdc_focus_times(tx_aper, 0, tx_delay(i,:));
       [rf_data_set{i}, t_start_set(i)] = calc_scat_multi(tx_aper, rx_aper, scat_pos, scat_amp);
   end
   [rfdata, t_start] = f_rf_comb(rf_data_set, t_start_set, fs);
   ```

5. **Free memory:**
   ```matlab
   xdc_free(tx_aper);
   xdc_free(rx_aper);
   ```

## Common Patterns

### Element Position Calculation

Matrix arrays use centered indexing:
```matlab
x_vec = (-(nx-1)/2 : (nx-1)/2) * pitch_x + center(1);
y_vec = (-(ny-1)/2 : (ny-1)/2) * pitch_y + center(2);
```

### Mathematical Element Subdivision

Standard pattern in all aperture functions:
```matlab
n_sub_x = ceil(physical_width / math_size);
sub_width = physical_width / n_sub_x;  % Actual subdivision size
```

### Physical Element Numbering

When concatenating apertures, physical element IDs must be sequential:
```matlab
rect = [phy_elem_id, corner1_xyz, corner2_xyz, corner3_xyz, corner4_xyz,
        weighting, width, height, center_xyz];
```

The first column contains physical element numbers (1-indexed in MATLAB, but Field-II uses 0-indexed internally).

## Language and Documentation

- Code contains mix of English and Chinese comments
- Function headers use English with MATLAB standard format
- Some newer functions (`xdc_rect_array.m`, `calc_rf.m`) use Chinese documentation
- Variable names are consistently English

## Known Limitations

- Rectangular elements must be aligned with x/y axes (no rotation support yet)
- No support for `xdc_triangles` interface (planned for arbitrary shapes)
- RCA array function (`xdc_rca_array.m`) is incomplete
- Multiple scattering ignored (Born approximation)
