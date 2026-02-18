#!/usr/bin/env bash

log()  { printf "[%s] INFO  %s\n"  "$(date '+%F %T')" "$*"; }
warn() { printf "[%s] WARN  %s\n"  "$(date '+%F %T')" "$*" >&2; }
die()  { printf "[%s] ERROR %s\n" "$(date '+%F %T')" "$*" >&2; exit 1; }