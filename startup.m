function startup
  cwd = fileparts(mfilename('fullpath'));

  add2path(fullfile(cwd, 'MAINS'), true);
  add2path(fullfile(cwd, 'MAGEKF'), true);
  add2path(fullfile(cwd, 'calibration'), true);
  add2path(fullfile(cwd, 'data'), true);
end

function add2path(p, doRecurse)
  if nargin<2 || ~doRecurse
    fprintf('Add path: %s\n', p);
    addpath(p);
  else
    fprintf('Recursively add path: %s\n', p);
    addpath(genpath(p));
  end
end
