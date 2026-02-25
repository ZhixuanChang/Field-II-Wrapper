# AGENTS.md

Guidelines for agentic coding agents working in this repository.

## Project Overview

Field-II-Wrapper is a MATLAB wrapper library for Field-II ultrasound simulation. It provides high-level functions for creating transducer apertures and managing simulation workflows. The library uses spatial impulse response methods to simulate ultrasound scattering with point scatterers under Born approximation.

## Build/Lint/Test Commands

This is a pure MATLAB project with no formal build system or package manager.

### Running Examples (Tests)
Examples serve as both documentation and tests. Run them in MATLAB:
```matlab
field_init  % Initialize Field-II (must be in MATLAB path)
example_linear_array_pwi      % Plane wave imaging with linear array
example_sparse_array_rect     % Sparse array 3D imaging
example_concatenate           % Aperture concatenation demo
example_rca_array_sector_scan % RCA array sector scanning
```

### Running a Single Test
There is no formal test framework. To test a function:
1. Open MATLAB
2. Ensure Field-II is in path and initialized: `field_init`
3. Run the example script or test file directly

### Linting
MATLAB has no standard linter. Use the MATLAB Editor's code analyzer (MLint) warnings as guidance.

## Code Style Guidelines

### File Naming
- Function files: `xdc_*.m` for aperture functions, `calc_*.m` for simulation functions, `f_*.m` for utility functions
- Example files: `example_*.m`
- Test files: `test_*.m`
- Template/temporary files: `tmp*.m`

### Function Header Format
Use MATLAB standard documentation format:
```matlab
function output = function_name(varargin)
% FUNCTION_NAME - Brief description (one line)
%
% Detailed description paragraph.
%
% Usage
%   output = function_name(param1, param2);
%
% Parameters
%   param1 - Description with units [unit].
%   param2 - Description.
%
% Return
%   output - Description.
%
% Example
%   result = function_name(1, 2);
%
% See also
%   related_function1, related_function2
%
% Author: Name (email)
% Date: Month DD, YYYY
% Version: X.X
% Comment: Version notes.
```

Chinese documentation is acceptable. Maintain consistency within each file.

### Variable Naming
- Use **snake_case** for function parameters and most variables: `elem_pos`, `num_elem`, `tx_delay`
- Use **camelCase** sparingly, mainly for legacy Field-II compatibility
- Constants and parameters at script top: lowercase with underscores
- Loop variables: `i`, `j`, `k` for simple loops; descriptive names for complex logic

### Code Structure
- Section markers: Use `%%` followed by section title for logical blocks
- Indentation: 4 spaces (MATLAB standard)
- Line length: Keep under 100 characters when practical

### Error Handling
Use MATLAB's `error()` function with descriptive messages:
```matlab
if nargin < 5
    error('function_name requires 5 input parameters');
end
if size(matrix, 2) ~= 3
    error('matrix must be N x 3');
end
```

### Optional Parameters
Use `varargin` with explicit argument count checking:
```matlab
if nargin == 8
    focus = [0, 0, 1000];
    center = [0, 0, 0];
elseif nargin == 10
    focus = varargin{9};
    center = [0, 0, 0];
else
    error('Invalid arguments number. Use help function_name.');
end
```

### Comments
- Comments are in **English or Chinese**
- Explain **why**, not **what** (code should be self-documenting)
- Document physical units in brackets: `[m]`, `[Hz]`, `[s]`
- Use block comments for complex algorithms

### Coordinate System
- X-axis: Element width direction
- Y-axis: Element height direction  
- Z-axis: Acoustic propagation direction (normal to transducer face)
- Arrays centered at origin by default

## Key Architectural Concepts

### Physical vs Mathematical Elements
Field-II subdivides physical elements into smaller "mathematical elements":
- Physical elements: Actual transducer elements (e.g., 64 elements)
- Mathematical elements: Computational subdivisions for accuracy
- `n_sub_x`, `n_sub_y` control subdivision density

### Aperture Functions Hierarchy
- **Low-level**: `xdc_rect_single.m`, `xdc_concatenate.m`
- **High-level**: `xdc_matrix_array.m`, `xdc_sparse_array_rect.m`, `xdc_rect_array.m`

### Important Constraints
- Rectangular elements must align with x/y axes (no rotation support yet)
- `xdc_rectangles()` requires: `rect` matrix, `center` matrix, `focus` vector
- Always call `xdc_free()` to release aperture memory after simulation

## Dependencies

- Field-II must be installed and in MATLAB path
- Initialize before use: `field_init`
- Key Field-II functions used: `xdc_rectangles`, `xdc_get`, `xdc_focus_times`, `xdc_apodization`, `calc_scat_multi`, `xdc_free`

## Memory Management

Always free apertures after simulation:
```matlab
xdc_free(tx_aper);
xdc_free(rx_aper);
```

When using `xdc_concatenate`, manually free the individual apertures in `aper_set` after concatenation.
