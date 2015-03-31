function [K, status] = PyMyAMI(varargin)

%   THIS IS A MODIFIED MATLAB/PYTHON INTERFACE FOR CALLING THE PITZER
%   COMPONENT OF THE MyAMI-V1.0 MODEL.
%
%   HOW TO CALL?
%   >>PITZERpath = '~/path/to/folder/of/MyAMI/code/Pitzer.py';
%   >>T = 20;
%   >>S = 35;
%   >>Ca = 0.0102821 * S/35;
%   >>Mg = 0.0528171* S/35;
%   >>B = xB*S / 1.80655 * 0.0000219;   
%   >>[K] = PyMyAMI(PITZERpath,num2str(T),num2str(S),num2str(Ca),num2str(Mg));
%
%   
%   Reference:
%   Hain,M.P., Sigman, D.M., Higgins, J.A., and Haug, G.H. (2015) 
%   The effects of secular calcium and magnesium concentration changes...
%   on the thermodynamics of seawater acid/base chemistry:...
%   Implications for Eocene and Cretaceous ocean carbon chemistry and buffering,
%   Global Biogeochemical Cycles, 29, doi:10.1002/2014GB004986


%   HEADER INFO FROM ORIGINAL MATLAB/PYTHON INTERFACE:
%   python Execute python command and return the result.
%   python(pythonFILE) calls python script specified by the file pythonFILE
%   using appropriate python executable.
%
%   python(pythonFILE,ARG1,ARG2,...) passes the arguments ARG1,ARG2,...
%   to the python script file pythonFILE, and calls it by using appropriate
%   python executable.
%
%   RESULT=python(...) outputs the result of attempted python call.  If the
%   exit status of python is not zero, an error will be returned.
%
%   [RESULT,STATUS] = python(...) outputs the result of the python call, and
%   also saves its exit status into variable STATUS.
%
%   If the python executable is not available, it can be downloaded from:
%     http://www.cpan.org
%
%   See also SYSTEM, JAVA, MEX.

%   Copyright 1990-2012 The MathWorks, Inc.
%   $Revision: 1.1.4.12 $
    
cmdString = '';

% Add input to arguments to operating system command to be executed.
% (If an argument refers to a file on the MATLAB path, use full file path.)
for i = 1:nargin
    thisArg = varargin{i};
    if ~ischar(thisArg)
        error(message('MATLAB:python:InputsMustBeStrings'));
    end
    if i==1
        if exist(thisArg, 'file')==2
            % This is a valid file on the MATLAB path
            if isempty(dir(thisArg))
                % Not complete file specification
                % - file is not in current directory
                % - OR filename specified without extension
                % ==> get full file path
                thisArg = which(thisArg);
            end
        else
            % First input argument is pythonFile - it must be a valid file
            error(message('MATLAB:python:FileNotFound', thisArg));
        end
    end

    % Wrap thisArg in double quotes if it contains spaces
    if isempty(thisArg) || any(thisArg == ' ')
        thisArg = ['"', thisArg, '"']; %#ok<AGROW>
    end

    % Add argument to command string
    cmdString = [cmdString, ' ', thisArg]; %#ok<AGROW>
end

% Execute python script
if isempty(cmdString)
    cmdString %screen output
    error(message('MATLAB:python:NopythonCommand'));
elseif ispc
    % PC
    pythonCmd = fullfile(matlabroot, 'sys\python\win32\bin\');
    cmdString = ['python' cmdString];
    pythonCmd = ['set PATH=',pythonCmd, ';%PATH%&' cmdString];
    [status, result] = dos(pythonCmd);
else
    % UNIX
    [status, ~] = unix('which python');
    if (status == 0)
        cmdString = ['python', cmdString];
        [status, result] = unix(cmdString);
    else
        error(message('MATLAB:python:NoExecutable'));
    end
end

% Check for errors in shell command
if nargout < 2 && status~=0
    error(message('MATLAB:perl:ExecutionError', result, cmdString)); % I rely on the PERL Error cataloge
end

K.Kw = zeros(3,1);
K.K1 = zeros(3,1);
K.K2 = zeros(3,1);
K.KsC = zeros(3,1);
K.Kb = zeros(3,1);
K.KsA = zeros(3,1);
K.K0 = zeros(3,1);
K.KSO4 = zeros(3,1);

K.Temp = 0;
K.Sal = 0;
K.Ca = 0;
K.Mg = 0;

[token, remain] = strtok(result);
% '[' parsed away

[token, remain] = strtok(remain);
K.Ca = str2num(token);
[token, remain] = strtok(remain);
K.Mg = str2num(token);
[token, remain] = strtok(remain);
K.Temp = str2num(token);
[token, remain] = strtok(remain);
K.Sal = str2num(token);

% empirical constants for modern composition Ca=0.0102821 & Mg=0.0528171
[token, remain] = strtok(remain);
K.Kw(1,1) = str2num(token);
[token, remain] = strtok(remain);
K.K1(1,1) = str2num(token);
[token, remain] = strtok(remain);
K.K2(1,1) = str2num(token);
[token, remain] = strtok(remain);
K.KsC(1,1) = str2num(token);
[token, remain] = strtok(remain);
K.Kb(1,1) = str2num(token);
[token, remain] = strtok(remain);
K.KsA(1,1) = str2num(token);
[token, remain] = strtok(remain);
K.K0(1,1) = str2num(token);
[token, remain] = strtok(remain);
K.KSO4(1,1) = str2num(token);

% Pitzer-model constants for [Ca] and [Mg] as specified
[token, remain] = strtok(remain);
K.Kw(2,1) = str2num(token);
[token, remain] = strtok(remain);
K.K1(2,1) = str2num(token);
[token, remain] = strtok(remain);
K.K2(2,1) = str2num(token);
[token, remain] = strtok(remain);
K.KsC(2,1) = str2num(token);
[token, remain] = strtok(remain);
K.Kb(2,1) = str2num(token);
[token, remain] = strtok(remain);
K.KsA(2,1) = str2num(token);
[token, remain] = strtok(remain);
K.K0(2,1) = str2num(token);
[token, remain] = strtok(remain);
K.KSO4(2,1) = str2num(token);

% Pitzer-model constants for modern composition Ca=0.0102821 & Mg=0.0528171
[token, remain] = strtok(remain);
K.Kw(3,1) = str2num(token);
[token, remain] = strtok(remain);
K.K1(3,1) = str2num(token);
[token, remain] = strtok(remain);
K.K2(3,1) = str2num(token);
[token, remain] = strtok(remain);
K.KsC(3,1) = str2num(token);
[token, remain] = strtok(remain);
K.Kb(3,1) = str2num(token);
[token, remain] = strtok(remain);
K.KsA(3,1) = str2num(token);
[token, remain] = strtok(remain);
K.K0(3,1) = str2num(token);
[token, remain] = strtok(remain);
K.KSO4(3,1) = str2num(token);



