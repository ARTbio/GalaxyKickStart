#!/usr/bin/env python

import argparse

from bioblend.galaxy import GalaxyInstance   

def create_user( url, email, password, username, key ):
    if email and password and username and key:
       galaxy_instance = GalaxyInstance(url, key)
       galaxy_instance.users.create_local_user(username,email,password)
    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="usage: python %prog [options]")

    parser.add_argument("-r", "--url",
                        default="http://localhost",
                        help="URL for the galaxy instance",)
    parser.add_argument("-u", "--username",
                        default="cloud",
                        help="Username to create. Defaults to 'cloud'",)
    parser.add_argument("-e", "--email",
                        default="cloud@galaxyproject.org",
                        help="Email for user",)
    parser.add_argument("-p", "--password",
                        default="password",
                        help="Password for user",)
    parser.add_argument("-k", "--key",
                        default="",
                        help="Master api key for the galaxy instance",)

    args = parser.parse_args()
    create_user(args.url, args.username, args.email, args.password, args.key)