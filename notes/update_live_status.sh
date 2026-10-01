#!/bin/bash
# Rewrites the "Live Cell Ranger status" block in notes/mpn-scrnaseq-project-status.md
# in place, between fixed markers, so it's safe to re-run repeatedly (idempotent).
set -euo pipefail
cd "$(dirname "$0")/.."
STATUS_FILE="notes/mpn-scrnaseq-project-status.md"
JOBID=18927392

running=$(squeue -u "$USER" -h -o "%i %T" | grep -c "^${JOBID}_.*RUNNING" || true)
pending=$(squeue -u "$USER" -h -o "%i %T" | grep -c "^${JOBID}_.*PENDING" || true)
done_count=$( (ls -d cellranger_count/*/outs 2>/dev/null || true) | wc -l | tr -d ' ')
failed=$(sacct -j "$JOBID" --format=JobID,State -n -P 2>/dev/null | grep -Ev '\.batch|\.extern' | grep -cE 'FAILED|CANCELLED|TIMEOUT|OUT_OF_ME' || true)
now=$(date '+%Y-%m-%d %H:%M:%S %Z')

completed_samples=""
if [ "$done_count" -gt 0 ]; then
  completed_samples=$( (ls -d cellranger_count/*/outs 2>/dev/null || true) | sed -E 's#cellranger_count/([^/]+)/outs#\1#' | paste -sd', ')
fi

BLOCK=$(cat <<EOF
<!-- LIVE-STATUS-START -->
## Live Cell Ranger status (auto-updated by notes/update_live_status.sh — safe to trust even if this session disconnects)
- Last checked: $now
- Job ${JOBID} (Stage 1, cellranger count): **${running} running, ${pending} pending, ${done_count}/19 completed, ${failed} failed**
- Completed samples: ${completed_samples:-none yet}
<!-- LIVE-STATUS-END -->
EOF
)

# Strip any prior block between markers, then append the fresh one.
if [ -f "$STATUS_FILE" ]; then
  awk '/<!-- LIVE-STATUS-START -->/{skip=1} !skip{print} /<!-- LIVE-STATUS-END -->/{skip=0}' "$STATUS_FILE" > "${STATUS_FILE}.tmp"
else
  : > "${STATUS_FILE}.tmp"
fi
printf '%s\n\n%s\n' "$(cat "${STATUS_FILE}.tmp")" "$BLOCK" > "$STATUS_FILE"
rm -f "${STATUS_FILE}.tmp"

echo "updated: running=$running pending=$pending done=$done_count failed=$failed at $now"
[ "$failed" -gt 0 ] && echo "FAILURE DETECTED: $failed failed task(s)"
if [ "$running" -eq 0 ] && [ "$pending" -eq 0 ]; then
  echo "ALL TASKS FINISHED (running=0 pending=0)"
  exit 0
fi
exit 42
