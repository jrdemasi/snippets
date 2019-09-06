#!/usr/bin/env python3

""" 
Used to trim ZFS snapshots on a FreeNAS appliance if they're being sent from another appliance
and rather than 1:1 replication you want a different retention policy on the backups.
There is nothing elegant about this, and in fact it's fairly inefficient.  It was a quick
solution to a problem I had.  PRs accepted. 
"""

import requests
import argparse

api_user = 'root'
api_pass = 'YourPassword'

def delete_snap(id):
    r = requests.delete("https://localhost/api/v1.0/storage/snapshot/{}".format(id), auth=(api_user, api_pass), verify=False)
    if r.status_code != 204:
        print("Could not delete snap: {}. HTTP response {}".format(id, r.status_code))
    return

parser = argparse.ArgumentParser()
parser.add_argument("--filesystem", help="The name of the filesystem to trim snaps for", required=True, type=str)
parser.add_argument("--keephourly", help="The number of hourly snaps to keep", required=True, type=int)
parser.add_argument("--keepdaily", help="The number of daily snaps to keep", required=True, type=int)
parser.add_argument("--keepweekly", help="The number of weekly snaps to keep", required=True, type=int)
parser.add_argument("--keepmonthly", help="The number of monthly snaps to keep", required=True, type=int)
args = parser.parse_args()

filesystem = args.filesystem
keep_hourly = args.keephourly
keep_daily = args.keepdaily
keep_weekly = args.keepweekly 
keep_monthly = args.keepmonthly

query_params = {"limit": 0}
r = requests.get("https://localhost/api/v1.0/storage/snapshot/", verify=False, auth=(api_user, api_pass), params=query_params)

all_snaps = r.json()

my_snaps = []
for snap in all_snaps:
    if filesystem == snap['filesystem']:
        my_snaps.append(snap['fullname'])

hourly_snaps = []
daily_snaps = []
weekly_snaps = []
monthly_snaps = []
for snap in my_snaps:
    if snap[-1:] == 'h':
        hourly_snaps.append(snap)
    elif snap[-1:] == 'd':
        daily_snaps.append(snap)
    elif snap[-1:] == 'w':
        weekly_snaps.append(snap)
    elif snap[-1:] == 'm':
        monthly_snaps.append(snap)
    else:
        print("No idea what's going on here")

# Sort all of our lists, useful later for trimming snaps to the proper numbers
hourly_snaps.sort()
daily_snaps.sort()
weekly_snaps.sort()
monthly_snaps.sort()

if len(hourly_snaps) > keep_hourly:
    to_delete = len(hourly_snaps) - keep_hourly
    snap_names = []
    counter = 0 
    while counter < to_delete:
        snap_names.append(hourly_snaps.pop(0))
        counter = counter + 1 
    for snap in snap_names:
        for result in all_snaps:
            if result['fullname'] == snap:
                delete_snap(result['id'])

if len(daily_snaps) > keep_daily:
    to_delete = len(daily_snaps) - keep_daily
    snap_names = []
    counter = 0 
    while counter < to_delete:
        snap_names.append(daily_snaps.pop(0))
        counter = counter + 1 
    for snap in snap_names:
        for result in all_snaps:
            if result['fullname'] == snap:
                delete_snap(result['id'])

if len(weekly_snaps) > keep_weekly:
    to_delete = len(weekly_snaps) - keep_weekly
    snap_names = []
    counter = 0 
    while counter < to_delete:
        snap_names.append(weekly_snaps.pop(0))
        counter = counter + 1 
    for snap in snap_names:
        for result in all_snaps:
            if result['fullname'] == snap:
                delete_snap(result['id'])

if len(monthly_snaps) > keep_monthly:
    to_delete = len(monthly_snaps) - keep_monthly
    snap_names = []
    counter = 0 
    while counter < to_delete:
        snap_names.append(monthly_snaps.pop(0))
        counter = counter + 1 
    for snap in snap_names:
        for result in all_snaps:
            if result['fullname'] == snap:
                delete_snap(result['id'])
