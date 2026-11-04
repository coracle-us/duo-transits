#!/usr/bin/env python3
"""Is the ASAS-SN Sky Patrol light-curve service up? Cheap check before running ground_check.py --asassn.

Prints the server's current message and block-server list.
Exit 0: block servers are listed (service looks up).
Exit 2: no block servers (light-curve downloads will fail; the service was down for maintenance from 2026 Oct 13).
Exit 3: the load balancer could not be reached.
"""
import sys

def main():
    from pyasassn.client import SkyPatrolClient
    try:
        c = SkyPatrolClient(verbose=False)
    except Exception as e:
        print(f"UNREACHABLE: {e}"); sys.exit(3)
    print(c.message.strip())
    print(f"block servers: {c.block_servers}")
    if not c.block_servers:
        print("DOWN"); sys.exit(2)
    print("UP")

if __name__ == "__main__":
    main()
