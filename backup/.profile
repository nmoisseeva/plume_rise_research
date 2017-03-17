cd 
# turn on colorized output
export CLICOLOR=1

# my prompt
export PS1='\[\033[0;32m\]\u@\h\[\033[0;36m\] \w\[\033[00m\]: '

# my aliases
alias python=python3.4
alias l='ls -CF'
alias ll='ls -AFhl'
alias tree='tree -F -I *~ -L 2'
alias wrfdir='cd ~/Applications/WRFV3/test/em_fire/'

alias df='df -h'
alias d='du -ch -d 1'
alias da='du -chsx *'
alias dl='du -d 1 | sort -rn | head -n 50'
alias ncl63='$NCARG_ROOT/ncarg6.3/bin/ncl'

#FDS Setting the environment for FDS and Smokeview.
# source ~/.bashrc_fds
alias smv='/Users/nadya2/Applications/FDS/FDS6/bin/smokeview'

#VAPOR
. /Applications/VAPOR/VAPOR.app/Contents/MacOS/vapor-setup.sh

# MacPorts path
export PATH=/Users/nadya2/Applications:/Users/nadya2/Applications/bin:/opt/local/bin:/opt/local/sbin:/opt/local/share/:/opt/local/lib:$PATH

# WRF/NCL paths
export FC=/opt/local/bin/gfortran-mp-5
export F90=/opt/local/bin/gfortran-mp-5
export CC=/opt/local/bin/gcc-mp-5
export NETCDF=/opt/local
export JASPERLIB=/opt/local/lib/
export JASPERINC=/opt/local/include/jasper/
export OMP_NUM_THREADS=7
# export BURF=1
# export CRTM=1
# export J='-j 4'
export NCARG_ROOT=/Users/nadya2/Applications
#export DAT_DIR=/Users/nadya2/WRFDA/testdata
export DYLD_FALLBACK_LIBRARY_PATH=/opt/local/lib/gcc5


# Added by Canopy installer on 2015-01-07
# VIRTUAL_ENV_DISABLE_PROMPT can be set to '' to make bashprompt show that Canopy is active, otherwise 1
#VIRTUAL_ENV_DISABLE_PROMPT=1 source /Users/nmoisseeva/Library/Enthought/Canopy_64bit/User/bin/activate

# virtualenvwrapper
# http://doughellmann.com/docs/virtualenvwrapper/tips.html
#export VIRTUALENVWRAPPER_PYTHON=/opt/local/bin/python
#export VIRTUALENVWRAPPER_VIRTUALENV=/opt/local/bin/virtualenv
#export WORKON_HOME=$HOME/.venvs
#export PROJECT_HOME=$HOME/Horizon/nowcasting
#source /opt/local/bin/virtualenvwrapper.sh

# virtualenv aliases
# http://blog.doughellmann.com/2010/01/virtualenvwrapper-tips-and-tricks.html
alias v='workon'
alias v.sw='workon'
alias v.de='deactivate'
alias v.mk='mkvirtualenv --no-site-packages'
alias v.mk_sitepackages='mkvirtualenv'
alias v.ad='add2virtualenv'
alias v.cd='cdvirtualenv'
alias v.cd_sitepackages='cdsitepackages'
alias v.ls='lssitepackages'
alias v.rm='rmvirtualenv'

# automatically activate a virtualenv when cd'ing into project directory
# http://justinlilly.com/python/virtualenv_wrapper_helper.html
has_virtualenv () {
    if [ -e .venv ]; then
        workon `cat .venv`
    fi
}
venv_cd () {
    builtin cd "$@" && has_virtualenv
}
alias cd='venv_cd'


