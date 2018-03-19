DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $DIR
export PATH=$DIR/bin:$PATH
export MROPATH=$DIR/mro
export GOBIN=$DIR/bin
export GOPATH=$DIR/lib/go
