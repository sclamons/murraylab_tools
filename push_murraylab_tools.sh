git add *
if [ "$1" != "" ]; then
    git commit -m "$1"
else
    git commit -m "Minor or otherwise undocumented automated reinstallation push"
fi
git push