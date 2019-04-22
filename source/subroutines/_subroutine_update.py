#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_update.py

# Import standard packages
import sys
import os
import subprocess

# Import included AddTag-specific modules
from . import subroutine

class ExtractParser(subroutine.Subroutine):
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'update'
        self.description = (
            "description:" "\n"
            "  Use the 'git' program to download the most recent version of" "\n"
            "  AddTag from the Git repository. Then update the software." "\n"
            "  May require your Atlassian login credentials." "\n"
        )
        self.help = "Update AddTag to the most recent version."
        self.epilog = (
            "example:" "\n"
            "  Try run AddTag with the following arguments:" "\n"
            "   $ python3 {__program__} {__subroutine__}" "\n"
        ).format(**dict(list(subroutine.__dict__.items()) + list({"__subroutine__": self.name}.items()))) # key:value pairs in the latter dict will replace any instances in the prior dict
        
        self.define_parser()
        self.define_arguments()
    
    def define_arguments(self):
        # Add mutually-exclusive arguments
        group = self.parser.add_mutually_exclusive_group()
        
        group.add_argument("--discard_local_changes", action="store_true",
            default=False, help="Force updated version to be a clean install.")
        
        group.add_argument("--keep_local_changes", action="store_true",
            default=False, help="Ensure the updated version does not remove local modifications.")
    
    def compute(self, args):
        """
        Perform the update
        """
        
        # Set the working directory as the AddTag folder.
        working_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))
        
        try:
            # Check to see it 'git' is available on PATH.
            command_list = ['git', 'rev-list', 'HEAD']
            cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
            #print(cp.stdout.decode().splitlines())
        except FileNotFoundError:
            # If 'git' is not on PATH, then exit with a non-zero status
            print("'git' program not found on PATH.", file=sys.stderr)
            sys.exit(1)
        
        # First, ask the user for confirmation (as their local files will may be re-written)
        user = input('Some local files may be removed or re-written. Proceed with update (y/n)? ')
        
        if (user in ['y', 'Y', 'yes', 'YES', 'Yes']):
            if args.discard_local_changes:
                # If you want the newest version, but you made local changes, then you can first discard your changes, and then update. Use the following two commands from inside the `addtag-project/` folder.
                command_list = ['git', 'reset', '--hard']
                cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
                command_list = ['git', 'pull']
                cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
            
            elif args.keep_local_changes:
                # Alternatively, if you want to keep the local modifications, you can use stash to hide them away before pulling, then reapply them afterwards.
                command_list = ['git', 'stash']
                cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
                command_list = ['git', 'pull']
                cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
                command_list = ['git', 'stash', 'pop']
                cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
            
            else:
                # If you would like to update your local copy to the newest version available, use the following command from within the `addtag-project/` directory.
                command_list = ['git', 'pull']
                cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
                
                # ./addtag update
                # Some local files may be removed or re-written. Proceed with update (y/n)? y
                # Password for 'https://tdseher@bitbucket.org': 
                # remote: Counting objects: 20, done.
                # remote: Compressing objects: 100% (20/20), done.
                # remote: Total 20 (delta 16), reused 0 (delta 0)
                # Unpacking objects: 100% (20/20), done.
                # From https://bitbucket.org/tdseher/addtag-project
                #    47ad5ff..920f9b0  master     -> origin/master
                # error: Your local changes to the following files would be overwritten by merge:
                #     source/donors.py
                #     source/nucleotides.py
                #     source/subroutines/_subroutine_confirm.py
                #     source/subroutines/_subroutine_generate.py
                # Please commit your changes or stash them before you merge.
                # Aborting
                # Traceback (most recent call last):
                #   File "./addtag", line 19, in <module>
                #     source.Main()
                #   File "/media/sf_VirtualBox_share/addtag-project/source/__init__.py", line 219, in __init__
                #     args.func(args)
                #   File "/media/sf_VirtualBox_share/addtag-project/source/subroutines/_subroutine_update.py", line 88, in compute
                #     cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
                #   File "/usr/lib/python3.5/subprocess.py", line 398, in run
                #     output=stdout, stderr=stderr)
                # subprocess.CalledProcessError: Command '['git', 'pull']' returned non-zero exit status 1

                
                
            print('Update finished.')
        else:
            print('No update performed.')
        
        # End 'compute()'
        