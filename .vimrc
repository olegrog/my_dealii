colorscheme ron
set hlsearch
set ts=2 sw=2 sts=2 et ai
filetype plugin on

"- Use the UTF-8 encoding
set encoding=utf-8
set fileencoding=utf-8

"- Return to last edit position when opening files
autocmd BufReadPost *
    \ if line("'\"") > 0 && line ("'\"") <= line("$") |
    \   exe "normal! g'\"" |
    \ endif

"- Highlight long lines
autocmd FileType python,c,cpp,sh set colorcolumn=80

"- Turn off autoindent
autocmd FileType * setlocal formatoptions-=c formatoptions-=r formatoptions-=o

"- Highlight trailing whitespaces
highlight ExtraWhitespace ctermbg=red guibg=red
match ExtraWhitespace /\s\+$/
autocmd BufWinEnter * match ExtraWhitespace /\s\+$/
autocmd InsertEnter * match ExtraWhitespace /\s\+\%#\@<!$/
autocmd InsertLeave * match ExtraWhitespace /\s\+$/
autocmd BufWinLeave * call clearmatches()

"- Highlight tabs
set list
set listchars=tab:▸\ ,

