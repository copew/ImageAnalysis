% silly function for committing to Git

function git_push_this(message)

!git add .
eval(['!git commit -m ''' message ''''])
!git push