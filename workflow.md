# Workflow Guide

This guide covers a variety of practices, mostly around Git, that I've picked up from software engineering. We're certainly not limited to these practices, but I've found they help avoid code conflict and keep things consistent.

## Some Basic Git Concepts

1. Repositories are top-level folders. They can have files or other folder inside of them.
2. Repositories exist on your local machine, but there is also a copy online in GitHub. This is like having a copy of a Word document on your computer while maintaining a copy in Google Drive. Your local version does not change the Drive version and vice versa unless you download or upload. We refer to the online version as the _remote_ or _origin_. Downloading is called _pulling_ or _fetching_ and uploading is called _pushing_. 
3. Saves in Git are called _commits_. Each commit is a snapshot of the entire repository at a given time.
4. Git relies on _branches_ for collaborative development. Each branch is a series of commits. The source of truth branch - the trunk of the tree, so to speak, is called `main`.

## Git Commands
These commands are usable from the command line, within the directory of a git repository.

- `git init` will initialize a git repository in a folder
- `git status` will show you what files you're tracking and what branch you're on, along with other goodies
- `git log` will show you your commit history
- `git add` will tell Git to track certain files
- `git add -p` will give you options to see all of the changes you're making before starting to track files
- `git add -A` will add all files that you don't have added in this folder
- `git commit -m "message here"` will save all of the files you're tracking as a snapshot in your commit history
- `git pull` pulls files from remote repositories into your local repository. Be careful to only do this on the `main` branch most of the time.
- `git push` pushes files from your local repository into remote repositories
- `git branch` shows you all branches and what branch you're currently on
- `git checkout -b new-branch-name` creates a new branch and switches to it
- `git checkout branch-name` switches to a given branch
- `git rebase branch-name` takes all the commits in your current branch and stacks them on top of the given branch's commit history

## A Guided Git Workflow
Suppose that you want to make some changes from your local repository.

Once you've got the repository cloned, `cd` into that directory. Run
```
git status
```
and you should see something like this:
```
On branch main
Your branch is up to date with 'origin/main'.

nothing to commit, working tree clean
```
If you have any local changes, they may also show up here, like this:
```
On branch main
Your branch is up to date with 'origin/main'.

Untracked files:
  (use "git add <file>..." to include in what will be committed)
	workflow.md

nothing added to commit but untracked files present (use "git add" to track)
```
To make edits, create a new branch like so:
```
git checkout -b my-branch-name
```
Make your edits and use the following commands to save as you work on your branch:
```
git add -A
git commit -m "Your commit message here"
```
Creating descriptive commit messages will help others to better understand your work.

Once you've committed all your edits, run the following series of commands:
```
git checkout main
git pull
git checkout my-branch-name
git rebase main
```
This will take your commits and save them on top of the remote main branch's commit history. Make sure that you're on your branch using `git checkout my-branch-name`. If this goes smoothly, run
```
git push origin my-branch-name
```
This pushes the local version of your branch up to GitHub, so that you'll have a remote version of your branch. Sometimes, however, this does not go as planned. You may have to fix merge conflicts.

To do so, use the command `git status`. This will show you which files have conflicts. You will then have to edit out the conflict messaging, choosing one version of the file to keep. After you have finished editing:
```
git add conflicted-file
```
Followed by:
```
git rebase --continue
```
And:
```
git push origin my-branch-name -f
```
Navigate to the repository on GitHub and you'll see a green pop-up that says "Create pull request."
See here for [more information](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests).

Create the pull request by writing a descriptive message and hitting the "Create Pull Request" button. At this point, you can send the pull request link to others to have them review your changes.

The pull request will also show you whether there are any conflicts with existing files in the `main` branch. Hopefully, these have already been fixed by your rebase! If there are conflicts, you can always start from the `git checkout main` up above and rebase once more.

Once you feel good about the code, hit the green merge button. This will merge your changes into main; the next time you run `git pull` on the `main` branch, you'll have the changes you made on your branch.


At this point, it's good practice to pull from main and delete your local branch:
```
git checkout main
git pull
git branch -d my-branch-name
```