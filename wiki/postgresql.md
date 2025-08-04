## ostgreSQL 서버 설치
- 이 가이드에서는 WSL2 Ubuntu에서의 설치 방법만 기술합니다.
- - 만약 Windows에 설치할려면 다음의 링크를 참고하시길 바랍니다. https://benkorea.github.io/NMIQ/posts/tools/PostgreSQL.html
- 또한 원래는 개인정보를 다루는 프로젝트와 공동으로 사용할 수 있도록 서비스계정으로 설치할려고 했지만 서비스계정의 uid와 gid를 컨테이너에 override 할 수 없어서 sudo 권한이 있는 개발자 계정으로 설치를 진행했습니다. - 사실 다음과 같은 검토를 했지만 현재까지 성공한 것이 없습니다.
- bitnami 이미지는 uid/gid 설정에 대한 환경변수가 없으므로 포기했습니다.
- sameersbn 이미지도 uid/gid 설정이 환경변수에 없었지만, 여기에서는 dockerfile로부터 build 하는 예시에 uid/gid를 사용자가 원하는대로 설정한 것이 있어 따라서 시도했지만 일단 저는 실패햇습니다. (github.com/sameersbn/docker-postgresql?tab=readme-ov-file)


## 전제조건
- Docker 및 Docker Compose가 설치되어 있어야 합니다
- 최소 2GB 메모리와 5GB 디스크 공간이 필요합니다
- 방화벽에서 포트 5432 열려있어야 합니다


## 설치 자동스트립트
- 프로젝트 루트에서 아래의 명령을 실행하시면 templates에 있던 docker-compose.yml괴 .env 파일이 /opt/docker/postgres로 복사되고, docker-compose up -d 명령이 실행됩니다.
```bash
sudo bash scripts/setup/wsl2_install_postgres.sh
```

## 컨테이너 내부의 id 확인
```bash
sudo docker exec postgres sh -c "id && id postgres"
```