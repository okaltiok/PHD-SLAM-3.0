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
    % Last Rev : 29/9/2022
    % Tested   : '9.8.0.1359463 (R2020a) Update 1'
    %
    % Copyright notice: You are free to modify, extend and distribute 
    %    this code granted that the author of the original code is 
    %    mentioned as the original author of the code.

    opt = 2; % datasets = {'Synthetic Data','Victoria Park'};
    
    clear mex
    
    % define paths and files
    folderParts = strsplit(pwd, filesep);
    switch opt
        case 1
            fprintf('Compiling MEX-files for synthetic data set\n')
            outdir = strjoin([folderParts(1:end-2) {'Synthetic Data' 'compiled_mex' ''}], filesep);
            model_path = strjoin({'-I.' 'Synthetic_Data_models' ''},filesep);
            model_file = strjoin({'Synthetic_Data_models' 'models.c'},filesep);
        case 2
            fprintf('Compiling MEX-files for Victoria Park data set\n')
            outdir = strjoin([folderParts(1:end-2) {'Victoria Park' 'compiled_mex' ''}], filesep);
            model_path = strjoin({'-I.' 'Victoria_Park_models' ''},filesep);
            model_file = strjoin({'Victoria_Park_models' 'models.c'},filesep);
    end

    linearAlgebra_path = strjoin({'-I.' ''},filesep);
    linearAlgebra_file = strjoin({'linearAlgebra.c'},filesep);
    murty_path = strjoin({'-I.' 'Shared_C++_Code' ''},filesep);
    murty_file = strjoin({'Shared_C++_Code' 'ShortestPathCPP.cpp'},filesep);
    

    % compile C/C++ files to Matlab MEX-files
    fprintf('compiling cost_matrix.c\n')
    compile_file = 'cost_matrix.c';
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir',outdir,...
        model_path, linearAlgebra_path, compile_file, model_file, linearAlgebra_file);
       
    fprintf('compiling gmiplid.c\n')
    compile_file = 'gmiplid.c';
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir',outdir,...
        model_path, linearAlgebra_path, compile_file, model_file, linearAlgebra_file);
    
    fprintf('compiling hypothesisReductionAlgorithm.c\n')
    compile_file = 'hypothesisReductionAlgorithm.c';
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir',outdir,...
        linearAlgebra_path, compile_file, linearAlgebra_file);
  
    fprintf('compiling kBest2DAssign.c\n')
    compile_file = 'kBest2DAssign.cpp';
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir',outdir,...
        murty_path, compile_file, murty_file);
    
    fprintf('compiling posterior_est.c\n')
    compile_file = 'posterior_est.c';
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir',outdir,compile_file);
    
    fprintf('compiling predict.c\n')
    compile_file = 'predict.c';
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir',outdir,...
        model_path, linearAlgebra_path, compile_file, model_file, linearAlgebra_file);
    
    fprintf('compiling sample.c\n')
    compile_file = 'sample.c';
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir',outdir,...
        linearAlgebra_path, compile_file, linearAlgebra_file);

    fprintf('compiling update.c\n')
    compile_file = 'update.c';
    mex('-R2018a','CXXFLAGS="$CXXFLAGS -std=c++17 -Wall"','-outdir',outdir,...
        model_path, linearAlgebra_path, compile_file, model_file, linearAlgebra_file);