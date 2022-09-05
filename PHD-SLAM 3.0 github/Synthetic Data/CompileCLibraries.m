function CompileCLibraries()

    % This function compiles the C/C++ files used with Synthetic data. If
    % you have compilation issues, read the following carefully!
    %
    % Check the compiler settings by calling:
    %
    % >>myCCompiler = mex.getCompilerConfigurations('C','Selected')
    % 
    % from your command window. The following settings have been used 
    % successfully:
    % 
    % myCCompiler = 
    % 
    %   CompilerConfiguration with properties:
    % 
    %              Name: 'MinGW64 Compiler (C)'
    %      Manufacturer: 'GNU'
    %          Language: 'C'
    %           Version: '6.3.0'
    %          Location: 'C:\ProgramData\MATLAB\SupportPackages\R2020a\3P.instrset\mingw_w64.instrset'
    %         ShortName: 'mingw64'
    %          Priority: 'E'
    %           Details: [1×1 mex.CompilerConfigurationDetails]
    %        LinkerName: 'C:\ProgramData\MATLAB\SupportPackages\R2020a\3P.instrset\mingw_w64.instrset\bin\gcc'
    %     LinkerVersion: ''
    %            MexOpt: 'C:\Apps\Matlab\R2020a-64bit\bin\win64\mexopts\mingw64.xml'
    %
    % If you don't have a compiler installed, go to the following website for
    % instructions on how to install it on your device.
    %
    % https://se.mathworks.com/help/matlab/matlab_external/install-mingw-support-package.html
    %
    % If you have multiple C or C++ compilers, choose MinGW by calling:
    %
    % >>mex -setup
    %
    % from your command window.
    % 
    %
    % Author   : Ossi Kaltiokallio
    %            Tampere University, Department of Electronics and
    %            Communications Engineering
    %            Korkeakoulunkatu 1, 33720 Tampere
    %            ossi.kaltiokallio@tuni.fi
    % Last Rev : 1/9/2022
    % Tested   : '9.8.0.1359463 (R2020a) Update 1'
    %
    % Copyright notice: You are free to modify, extend and distribute 
    %    this code granted that the author of the original code is 
    %    mentioned as the original author of the code.

    clear mex

    fprintf('compiling cost_matrix.c\n')
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir','./compiled_mex',...
        '-I../Shared Files/mex source/Synthetic_Data_models/','../Shared Files/mex source/cost_matrix.c',...
        '../Shared Files/mex source/Synthetic_Data_models/models.c');

    fprintf('compiling hypothesisReductionAlgorithm.c\n')
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir','./compiled_mex',...
        '-I../Shared Files/mex source/Synthetic_Data_models/','../Shared Files/mex source/hypothesisReductionAlgorithm.c',...
        '../Shared Files/mex source/Synthetic_Data_models/models.c');

    fprintf('compiling kBest2DAssign.c\n')
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir','./compiled_mex',...
        '-I../Shared Files/mex source/Shared_C++_Code/','../Shared Files/mex source/kBest2DAssign.cpp',...
        '../Shared Files/mex source/Shared_C++_Code/ShortestPathCPP.cpp');

    fprintf('compiling mhrbiploid.c\n')
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir','./compiled_mex',...
        '-I../Shared Files/mex source/Synthetic_Data_models/','../Shared Files/mex source/mhrbiploid.c',...
        '../Shared Files/mex source/Synthetic_Data_models/models.c');

    fprintf('compiling posterior_est.c\n')
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir','./compiled_mex',...
        '../Shared Files/mex source/posterior_est.c');

    fprintf('compiling predict.c\n')
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir','./compiled_mex',...
        '-I../Shared Files/mex source/Synthetic_Data_models/','../Shared Files/mex source/predict.c',...
        '../Shared Files/mex source/Synthetic_Data_models/models.c');

    fprintf('compiling sample.c\n')
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir','./compiled_mex',...
        '-I../Shared Files/mex source/Synthetic_Data_models/','../Shared Files/mex source/sample.c',...
        '../Shared Files/mex source/Synthetic_Data_models/models.c');

    fprintf('compiling update.c\n')
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir','./compiled_mex',...
        '-I../Shared Files/mex source/Synthetic_Data_models/','../Shared Files/mex source/update.c',...
        '../Shared Files/mex source/Synthetic_Data_models/models.c');
