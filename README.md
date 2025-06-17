# SAR Algorithm

This repository was migrated from Gitee to GitHub. It contains the implementation of the SAR (Synthetic Aperture Radar) algorithm.

## Introduction

This project is a course assignment for SAR signal processing and motion compensation. It implements several key SAR imaging algorithms including Range-Doppler Algorithm (RDA), Chirp Scaling Algorithm (CSA), and wavenumber domain algorithm (wKA).

## Features

- Implementation of multiple SAR imaging algorithms:
  - Range-Doppler Algorithm (RDA)
  - Chirp Scaling Algorithm (CSA)
  - Wavenumber domain algorithm (wKA)
- Image enhancement capabilities
- Point target analysis
- Performance metrics calculation (IRW, PSLR, ISLR)
- Lee filter implementation for speckle reduction

## Project Structure

- :open_file_folder:src
  - :open_file_folder:functions
    - :page_facing_up:RDA.m - Range-Doppler Algorithm implementation
    - :page_facing_up:CSA.m - Chirp Scaling Algorithm implementation
    - :page_facing_up:wKA.m - Wavenumber domain algorithm implementation
    - :page_facing_up:read_data.m - Data reading utilities
    - :page_facing_up:generate_point_data.m - Point target data generation
    - :page_facing_up:f_point_analyze.m - Point target analysis
    - :page_facing_up:f_IRW_PSLR_ISLR.m - Performance metrics calculation
    - :page_facing_up:Lee_filter.m - Speckle reduction filter
    - ...
  - :page_facing_up:main.m - Main execution file
  - :page_facing_up:RDA_sim.m - RD algorithm simulation script
  - :page_facing_up:CS_sim.m - RD algorithm simulation script
  - :page_facing_up:wKA_sim.m - RD algorithm simulation script
  - :page_facing_up:image_enhance.m - Image enhancement script
- :open_file_folder:test - Test files and examples
- :open_file_folder:trash - Temporary files
- :page_facing_up:变量名约定.md - Variable naming conventions

## Requirements

- MATLAB (R2018b or later recommended)
- Signal Processing Toolbox
- Image Processing Toolbox

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/sar-algorithm.git
```

2. Add the project directory to your MATLAB path:
```matlab
addpath(genpath('path/to/sar-algorithm'));
```

## Usage

This repository contains implementations of three SAR imaging algorithms: RDA, CSA, and wKA. To use these algorithms:

1. Add the project directory to your MATLAB path
2. Run the corresponding simulation script:
```matlab
RDA_sim  % Range-Doppler Algorithm
CS_sim   % Chirp Scaling Algorithm
wKA_sim  % Wavenumber domain algorithm
```

For detailed parameter settings and examples, please refer to the comments in each script file.

## Contributing

We welcome contributions! If you find any issues or have suggestions for improvements:

1. Fork the repository
2. Create your feature branch
3. Commit your changes
4. Push to the branch
5. Create a Pull Request

## Authors

### Original Authors (Gitee)
This project was originally developed and maintained by:
- Hongxiang Li@[lihongxiang1998](https://gitee.com/lihongxiang1998)
- Jianwei Liu@[XiaoLiu2021](https://gitee.com/XiaoLiu2021)
- Fei Zhao@[zhaofei2048](https://gitee.com/zhaofei2048)

### Current Maintainers
This repository is currently maintained by me after migration from Gitee to GitHub.

## License

This project is open source and available for any use. While not required, we appreciate if you reference this repository when using the code.

## Acknowledgments

This repository was migrated from Gitee to GitHub. The maintainers have made minor code corrections and formatting adjustments during the migration process. When citing or acknowledging this work, please refer to the original authors listed above.

