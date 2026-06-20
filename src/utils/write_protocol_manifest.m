function manifest_file = write_protocol_manifest(out_dir, protocol, seed_resolution)
% WRITE_PROTOCOL_MANIFEST  Write protocol and seed lineage as JSON.

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
payload = struct();
payload.Protocol = protocol;
payload.SeedResolution = seed_resolution;
payload.WrittenAt = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

manifest_file = fullfile(out_dir, 'protocol_manifest.json');
fid = fopen(manifest_file, 'w');
if fid < 0
    error('write_protocol_manifest:openFailed', ...
        'Unable to write protocol manifest: %s', manifest_file);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', jsonencode(payload, 'PrettyPrint', true));
end

