#!/bin/bash

set -e

# Docker 권한 확인
docker ps >/dev/null || { echo "Docker 권한 없음. 'sudo usermod -aG docker $USER' 실행 후 재로그인"; exit 1; }

# 설치
sudo cp -r "$(dirname "$0")/../../templates/postgres" /opt/
sudo mkdir -p /opt/postgres/data
sudo chown -R $USER:$USER /opt/postgres

# 실행
cd /opt/postgres && docker compose up -d

echo "완료: localhost:5432, DB:ai4ref, User:postgres, Pass:postgres"