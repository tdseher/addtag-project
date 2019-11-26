#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_update.py

# Import standard packages
import sys
import os
import shutil
import subprocess

# Import included AddTag-specific modules
from . import subroutine
from .subroutine import __repository__

class UpdateParser(subroutine.Subroutine):
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'update'
        self.description = (
            "description:" "\n"
            "  Use the 'git' program to download the most recent version of" "\n"
            "  AddTag from the Git repository. Then update the software." "\n"
            "  If the Git repository is private, it may request your login credentials." "\n"
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
            cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except FileNotFoundError:
            # If 'git' is not on PATH, then exit with a non-zero status
            print("ERROR: 'git' program not found on PATH. No update performed.", file=sys.stderr)
            sys.exit(1)
        except subprocess.CalledProcessError as error:
            #print('ERROR" {}'.format(error.stderr.decode()))
            # This error occurrs if '.git' folder does not exist.
            # If it doesn't, then the source should be checked out
            print("WARNING: '.git' folder does not exist.')")
            
            command_list = ['git', 'clone', __repository__, 'update']
            print("Running: '{}'".format(' '.join(command_list)))
            cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
            
            print("Moving '.git' folder")
            shutil.move('update/.git', '.git')
            print("Removing 'update' folder")
            shutil.rmtree('update')
        
        # Obtain current version information
        command_list = ['git', 'log', '-n', '1', 'HEAD', '--format=%H,%cd']
        print("Running: '{}'".format(' '.join(command_list)))
        current_hash, current_date = None, None
        try:
            cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            current_hash, current_date = cp.stdout.decode().splitlines()[0].split(',', 1)
        except subprocess.CalledProcessError as error:
            if error.stderr:
                print('ERROR: {}'.format(error.stderr.decode()))
        
        # Run 'git fetch' to "safely" download the repo from the remote repository
        command_list = ['git', 'fetch']
        print("Running: '{}'".format(' '.join(command_list)))
        try:
            cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as error:
            if error.stderr:
                print('ERROR: {}'.format(error.stderr.decode()))
        
        # Obtain updated version information
        updated_hash, updated_date = None, None
        command_list = ['git', 'log', '-n', '1', 'FETCH_HEAD', '--format=%H,%cd']
        print("Running: '{}'".format(' '.join(command_list)))
        try:
            cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            updated_hash, updated_date = cp.stdout.decode().splitlines()[0].split(',', 1)
        except subprocess.CalledProcessError as error:
            if error.stderr:
                print('ERROR: {}'.format(error.stderr.decode()))
        
        if (current_hash and updated_hash and (current_hash != updated_hash)):
            print('Current version: {} ({})'.format(current_hash[:7], current_date))
            print('New version: {} ({})'.format(updated_hash[:7], updated_date))
            print('')
            
            # Compare current version with remote/fetched version
            command_list = ['git', 'log', '--oneline', '--decorate=no', 'HEAD..FETCH_HEAD']
            print("Running: '{}'".format(' '.join(command_list)))
            cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
            changes = cp.stdout.decode().splitlines()
            print('Change log:')
            for change in changes:
                print(' * {}'.format(change.replace('[skip ci]', '')))
            print('')
            
            # First, ask the user for confirmation (as their local files will may be re-written)
            user = input('Some local files may be removed or re-written. Proceed with update (y/n)? ')
            
            if (user in ['y', 'Y', 'yes', 'YES', 'Yes']):
                if args.discard_local_changes:
                    # If you want the newest version, but you made local changes, then you can first discard your changes,
                    # and then update. Use the following two commands from inside the `addtag-project/` folder.
                    #command_list = ['git', 'reset', '--hard', 'HEAD']
                    #print("Running: '{}'".format(' '.join(command_list)))
                    #cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
                    
                    # Hide away the modifications
                    command_list = ['git', 'stash', '--include-untracked']
                    print("Running: '{}'".format(' '.join(command_list)))
                    cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
                    
                    # Remove the hidden modifications
                    command_list = ['git', 'stash', 'drop']
                    print("Running: '{}'".format(' '.join(command_list)))
                    cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
                    
                    # Download the new version
                    command_list = ['git', 'pull']
                    print("Running: '{}'".format(' '.join(command_list)))
                    cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
                
                elif args.keep_local_changes:
                    # Alternatively, if you want to keep the local modifications, you can use stash to hide them away
                    # before pulling, then reapply them afterwards.
                    command_list = ['git', 'stash', '--include-untracked']
                    print("Running: '{}'".format(' '.join(command_list)))
                    cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
                    
                    # Download the new version
                    command_list = ['git', 'pull']
                    print("Running: '{}'".format(' '.join(command_list)))
                    cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
                    
                    # Restore hidden modifications
                    command_list = ['git', 'stash', 'pop']
                    print("Running: '{}'".format(' '.join(command_list)))
                    try:
                        cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    except subprocess.CalledProcessError as error:
                        if error.stderr:
                            print('ERROR: {}'.format(error.stderr.decode()))
                            print('Conflicts exist between the new version and the local modifications.')
                            print("Please inspect the 'stash' to manually resolve these conflicts.")
                            print("Notwithstanding, you have newest version of AddTag.")
                    
                    # Consider adding the following approach in case there are conflicts:
                    #   git stash show -p | git apply
                    #   git stash drop
                    # or
                    #   git stash show -p | git apply --3way
                    #   git stash drop
                
                else:
                    # If you would like to update your local copy to the newest version available, use the following
                    # command from within the `addtag-project/` directory.
                    try:
                        command_list = ['git', 'pull']
                        print("Running: '{}'".format(' '.join(command_list)))
                        cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    except subprocess.CalledProcessError as error:
                        if error.stderr:
                            print('ERROR: {}'.format(error.stderr.decode()))
                        print("ERROR: 'git pull' returned non-zero exit status.", file=sys.stderr)
                        print("You may have local changes that must either be discarded or kept.", file=sys.stderr)
                        print("Please see 'addtag update --help' for more information.", file=sys.stderr)
                        sys.exit(1)
                    
                print('Update finished.')
            else:
                print('No update performed.')
        
        else:
            print('No updated version available.')
        
        # End 'compute()'
