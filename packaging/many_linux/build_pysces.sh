#!/bin/bash

adduser -u 1000 jr
su -c /io/packaging/many_linux/build_script.sh jr
