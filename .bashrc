# Some of bash colors
export RED='\033[0;31m'
export CYAN='\033[0;36m'
export BBLUE='\033[1;34m'
export NC='\033[0m'

# Command line
_parse_error_code() {
    local code=$?
    [[ $code -gt 0 ]] || return
    echo "[$code] "
}
# \[..color..\] should be used to preserve correct shell line width
export PS1="\[$RED\]\$(_parse_error_code)\[$CYAN\]\h\[$NC\]:\[$BBLUE\]\w\[$NC\]\$ "

# Bash history
export HISTTIMEFORMAT="%d/%m/%y %T "
export HISTSIZE=10000
export HISTFILESIZE=10000
export HISTIGNORE="ls:ll:cd:cd -"
export HISTCONTROL="ignorespace:erasedups"
shopt -s histappend

# BUG: history is not loaded by default
export HISTFILE=$HOME/.bash_history
history -r

# Colorize less output
export LESSOPEN="| $(which source-highlight) -i %s -f esc256 --style-file esc256.style"
export LESS=" -R"

# Find files in the deal.II distribution and among the local source code
_find_dealii_files() {
    local res1 res2
    res1="$(locate "/$1" | grep "/usr/src/dealii" | grep -v tests)"
    res2="$(find ~ -wholename "*/$1*" | grep -v examples)"
    if [[ -z "$res1$res2" ]]; then
        echo -e "${RED}Missing filename"\!"${NC}" >&2
    else
        echo -e "$res1\n$res2"
    fi
}

# Open an deal.II file
s() {
    local files
    files="$(_find_dealii_files "$1")"
    [[ -z "$files" ]] && return 1
    if [[ "$(wc -w <<< "$files")" -gt 1 ]]; then
        echo -e "${RED}There are several files with the same name." >&2
        echo -e "Specify the unique path using the parent directory.${NC}" >&2
        echo "$files" | tr ' ' '\n' >&2
        return 2
    fi
    local file=${files/$'\n'}
    if [[ -w "$file" ]]; then
        vim "$file"
    else
        less -R "$file"
    fi
}

# List files in an deal.II directory
d() {
    local files
    declare -A dirs
    files="$(_find_dealii_files "$1")"
    for file in $files; do
        if [[ -f "$file" && "$file" =~ "$1"$ ]]; then
            dirs["$(dirname "$file")"]=1
        elif [[ -d "$file" && "$file" =~ "$1"$ ]]; then
            dirs["$file"]=1
        fi
    done
    for dir in "${!dirs[@]}"; do
        echo -e "Directory: ${BBLUE}$dir${NC}"
        ls "$dir"
    done
}

# Aliases
alias grep='grep --color=auto'
alias ls='ls --color=auto -G'
alias ll='ls -lFh'
alias diff='colordiff'
