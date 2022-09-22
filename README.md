My deal.II solvers and problems
========

The latest version of [deal.II](https://www.dealii.org) is supported.

Installation
-------
The repository can be cloned by
```
    git clone https://github.com/olegrog/dealii.git
```
To build an augmented Docker image, run
```
    cd dealii
    docker build -f Dockerfile -t my_dealii .
```
The new image can be used within an intermediate container as
```
    docker run -it --rm -u="$(id -u):$(id -g)" -v="$(pwd):/home/dealii" my_dealii
```
