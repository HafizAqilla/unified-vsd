function reset_matlab_startup_path()
% RESET_MATLAB_STARTUP_PATH
% -----------------------------------------------------------------------
% Resets the saved MATLAB startup path to a clean user-level pathdef that
% does not recursively pin this repository or transient worktrees.
%
% This avoids startup warnings from stale repo-relative folders such as
% "lib/" and "modules/" that do not exist at the project root.
%
% INPUTS:
%   None.
%
% OUTPUTS:
%   None. Writes a user-level pathdef.m under userpath().
%
% ASSUMPTIONS:
%   - Project runtime adds its own clean path in main_run.m.
%   - Keeping repository subfolders out of MATLAB startup path is safe.
%
% REFERENCES:
%   [1] MATLAB savepath / userpath documentation.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-30
% VERSION:  1.0
% -----------------------------------------------------------------------

root_dir = fileparts(fileparts(mfilename('fullpath')));
user_startup_dir = strtrim(userpath);

if isempty(user_startup_dir)
    error('reset_matlab_startup_path:noUserPath', ...
        'MATLAB userpath() is empty; cannot write user-level pathdef.');
end

if ~isfolder(user_startup_dir)
    mkdir(user_startup_dir);
end

pathdef_target = fullfile(user_startup_dir, 'pathdef.m');
backup_target = fullfile(user_startup_dir, ...
    sprintf('pathdef_backup_%s.m', datestr(now, 'yyyymmdd_HHMMSS')));

if exist(pathdef_target, 'file') == 2
    copyfile(pathdef_target, backup_target);
    fprintf('[reset_matlab_startup_path] Backed up existing user pathdef to:\n  %s\n', ...
        backup_target);
end

current_paths = regexp(path, pathsep, 'split');
current_paths = current_paths(~cellfun('isempty', current_paths));

repo_prefix = [regexptranslate('escape', root_dir) '($|[\\/])'];
keep_mask = cellfun(@(p) isempty(regexp(p, repo_prefix, 'once')), current_paths);
clean_paths = current_paths(keep_mask);

path(strjoin(clean_paths, pathsep));
save_status = savepath(pathdef_target);

if save_status ~= 0
    error('reset_matlab_startup_path:saveFailed', ...
        'savepath failed for user-level pathdef: %s', pathdef_target);
end

fprintf('[reset_matlab_startup_path] Saved clean startup path to:\n  %s\n', ...
    pathdef_target);
fprintf('[reset_matlab_startup_path] Removed %d repo-specific startup entries.\n', ...
    sum(~keep_mask));

end
