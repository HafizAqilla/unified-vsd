# CHECK_NOISE_PROGRESS
# -----------------------------------------------------------------------
# Progress of a measurement-noise batch (run_noise_simulation.m), readable
# while MATLAB is still running. Read-only: it never touches the batch.
#
# USAGE (from the repository root):
#   .\scripts\check_noise_progress.ps1
#   .\scripts\check_noise_progress.ps1 -Workers 4 -Watch      # refresh every 5 min
#
# A replicate is DONE when replicates\mcNNN.mat exists (completed or failed;
# the status itself is in summary\noise_runs.csv after aggregation). It is
# RUNNING when its run folder exists but the .mat does not yet.
# -----------------------------------------------------------------------
param(
    [string]$BatchDir = "results\batches\noise_20260927_134251",
    [int]$Workers = 3,
    [double]$MinutesPerReplicate = 84,
    [switch]$Watch
)

function Show-Progress {
    $manifest = Get-Content (Join-Path $BatchDir "manifest.json") -Raw | ConvertFrom-Json
    $M = [int]$manifest.num_replicates
    $repDir = Join-Path $BatchDir "replicates"

    $done = @(); $running = @(); $pending = @()
    for ($k = 1; $k -le $M; $k++) {
        $tag = "mc{0:D3}" -f $k
        if (Test-Path (Join-Path $repDir "$tag.mat")) { $done += $k }
        elseif (Get-ChildItem $repDir -Directory -Filter "*_$tag" -ErrorAction SilentlyContinue) { $running += $k }
        else { $pending += $k }
    }

    $pct = [math]::Round(100 * $done.Count / $M)
    $bar = ('#' * [math]::Floor($pct / 5)).PadRight(20, '.')
    $left = $running.Count + $pending.Count
    $etaHours = [math]::Round([math]::Ceiling($left / $Workers) * $MinutesPerReplicate / 60, 1)

    Clear-Host
    Write-Host ("Noise batch  {0}" -f $BatchDir)
    Write-Host ("Checked      {0}" -f (Get-Date -Format "yyyy-MM-dd HH:mm"))
    Write-Host ""
    Write-Host ("[{0}] {1}%   {2} of {3} done" -f $bar, $pct, $done.Count, $M)
    Write-Host ("Done     : {0}" -f ($(if ($done) { $done -join ', ' } else { '-' })))
    Write-Host ("Running  : {0}" -f ($(if ($running) { $running -join ', ' } else { '-' })))
    Write-Host ("Pending  : {0}" -f $pending.Count)
    Write-Host ("Estimate : about {0} h left ({1} workers, ~{2} min per replicate)" -f $etaHours, $Workers, $MinutesPerReplicate)

    $last = Get-ChildItem $repDir -Filter "mc*.mat" -ErrorAction SilentlyContinue |
        Sort-Object LastWriteTime | Select-Object -Last 1
    if ($last) { Write-Host ("Last finished: {0} at {1}" -f $last.BaseName, $last.LastWriteTime.ToString("HH:mm")) }

    $log = Join-Path $BatchDir "run_noise_simulation.log"
    if (Test-Path $log) {
        Write-Host ""
        Write-Host "Last log lines:"
        Get-Content $log -Tail 5 | ForEach-Object { Write-Host ("  " + $_) }
    }
}

if ($Watch) {
    while ($true) { Show-Progress; Start-Sleep -Seconds 300 }
} else {
    Show-Progress
}
