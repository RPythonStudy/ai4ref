#!/bin/bash

set -e

# Docker 권한 확인
docker ps >/dev/null || { echo "Docker 권한 없음. 'sudo usermod -aG docker $USER' 실행 후 재로그인"; exit 1; }

# 설치
sudo mkdir -p /opt/docker
sudo cp -r "$(dirname "$0")/../../templates/postgres" /opt/docker/
sudo mkdir -p /opt/docker/postgres/data
sudo chown -R 1000:1000 /opt/docker/postgres

# 실행
cd /opt/docker/postgres && docker compose up -d

echo "완료: localhost:5432, DB:ai4ref, User:postgres, Pass:postgres"