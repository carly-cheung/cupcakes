# Cupcakes and coding 2017-10-12

Today we did:

- Installed `zsh` for the terminal: http://ohmyz.sh/
  - if it wasn't working, we fixed it by activating the editor using this command: `chsh -s /bin/zsh`
- Installed atom: https://atom.io/
  - Added the command line: http://flight-manual.atom.io/getting-started/sections/atom-basics/
- Installed homebrew: https://brew.sh/
  - `brew install tree` to see the folder structure

## Answers to absolute vs relative paths quiz

http://swcarpentry.github.io/shell-novice/02-filedir/

1. `cd .` keeps you in the same place
2. `cd /` takes you to the "root" folder
  - "root" is like the beginning of time for the computer
3. `cd /home/amanda`
4. `cd ../..` takes you up two folders from your current location. For Amanda, it takes her back to `/Users` which is not her home directory
5. `cd ~` Takes you back to home.
6. `cd home` gives you an error message because there's no such directory called "home"
1. `cd ~/data/..` Takes first to home, then data, a level back up, so yes it takes her back home.
1. `cd` Takes you back home. Also absolute (doesn't matter where you are)
1. `cd ..` Takes her up one hierarchical folder up so she'll be in `/Users/amanda`
