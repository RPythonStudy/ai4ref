# PostgreSQL 서버 설치
- 이 프로젝트 자체는 개인정보를 다루지 않으므로 PostgreSQL 서버를 서비스계정으로 설치할 필요는 없습니다. 그러나 PostgreSQL은 개인정보를 다루는 AI4RM 등의 프로젝트와 함께 사용할 수 있습니다. 이러한 측면에서는 개인정보를 다루는 프로젝트와 공동으로 사용할 수 있도록 서비스계정으로 설치를 합니다.
- PostgreSQL 회사에서 제공하는 공식 Docker 이미지는 컨테이너 내부 운영체제의 계정 uid/gid를 조정할 수 없습니다. 따라서 bitnami 도커 이미지를 사용합니다. 

## 전제조건
- Docker 및 Docker Compose가 설치되어 있어야 합니다
- 최소 2GB 메모리와 5GB 디스크 공간이 필요합니다
- 방화벽에서 포트 5432 열려있어야 합니다

## postgres 계정 생성
- 보안을 강화하기 위해 서비스 전용 계정을 생성하는 것을 권장합니다.
```bash
sudo useradd postgres
```

## sudoers 설정
- 일반적으로 서비스계정을 docker group에 포함을 해도 보안에 바람직합니다.
- 그러나 보안을 더욱 강화하기 위해 sudoers를 사용할 수도 있습니다.
```bash
sudo visudo -f /etc/sudoers.d/postgres-docker
```
- 아래를 추가하여 postgres가 sudo 명령을 사용할 수 있도록 합니다.
```plaintext
postgres ALL=(ALL) NOPASSWD: /usr/bin/docker
```


## 디렉토리 생성
- PostgreSQL 설치 디렉토리를 생성합니다.
```bash
sudo mkdir -p /opt/docker/postgres
```

## 디렉토리 소유권 설정
```bash
sudo chown postgres:postgres /opt/docker/postgres
```

## 디렉토리 권한 설정
- 프로젝트에서는 breakglassadmin 계정을 고려해서 750 권한을 사용합니다.
```bash
sudo chmod -R 750 /opt/docker/postgres
```

## 계정변경
- bitwarden 계정에서 설치를 진행하기 위해 계정을 전환합니다.
```bash 
su - bitwarden
```

## 설치 대상 폴더로 이동
```bash
cd /opt/ai4rm/docker/bitwarden
```

## 설치 스크립트 다운로드
```bash
curl -Lso bitwarden.sh "https://func.bitwarden.com/api/dl/?app=self-host&platform=linux" && chmod 700 bitwarden.sh
```

## sudo로 설치 진행
- 개발사 매뉴얼은 ./bitwarden.sh install 이지만, AI4RM 프로젝트에서는 아래와 같이 sudo로 설치합니다.
```bash
sudo ./bitwarden.sh install
```
## 설치 중 옵션
- 개발사의 매뉴얼에 따라 진행하시면 됩니다. (https://bitwarden.com/help/install-on-premise-linux/)

## user override 설정
- bitwarden.sh 설치 후, /opt/ai4rm/docker/bitwarden/bwdata/docker/docker-compose.override.yml 파일을 생성하여 bitwarden 컨테이너의 사용자와 그룹을 설정합니다.
- 아래의 1001:1001은 예시일 뿐이며, 위에서 id 명령어로 확인한 host의 uid:gid를 사용해야 합니다.
```yaml
# 파일 경로: /opt/ai4rm/docker/bitwarden/bwdata/docker/docker-compose.override.yml
services:
  api:
    user: "1001:1001"
  web:
    user: "1001:1001"
  admin:
    user: "1001:1001"
  identity:
    user: "1001:1001"
  sso:
    user: "1001:1001"
  events:
    user: "1001:1001"
  notifications:
    user: "1001:1001"
  attachments:
    user: "1001:1001"

```
- 이 파일은 프로젝트 폴더 templates/bitwarden/docker-compose.override.yml에서 복사하여 사용할 수 있습니다. 
- docker-compose.override.yml 복사는 권한문제 때문에 bitwarden 계정에서 진행합니다.
```bash
cp /opt/ai4rm/templates/bitwarden/docker-compose.override.yml /opt/ai4rm/docker/bitwarden/bwdata/docker/docker-compose.override.yml
```
- 소유계정설정 진행해야 합니다.
```bash
sudo chown bitwarden:bitwarden /opt/ai4rm/docker/bitwarden/bwdata/docker/docker-compose.override.yml
```
- 권한설정
```bash
sudo chmod 644 /opt/ai4rm/docker/bitwarden/bwdata/docker/docker-compose.override.yml
```


# CLI 설치
- CLI는 Bitwarden 서버와 상호작용하기 위해 프로젝트 내에서 호출하므로 아래 단계를 따라 설치합니다.

## 시스템 전역 설치
```bash
sudo curl -L https://vault.bitwarden.com/download/?app=cli&platform=linux -o /usr/local/bin/bw
```

## 실행 권한 부여
```bash
sudo chmod +x /usr/local/bin/bw
```

## 설치 확인
```bash
bw --version
```

## 보안 고려사항 및 감사 로깅
- 이부분은 개념상의 예시이며 수정이 필요합니다.

### Bitwarden 설치 감사 로그
- AI4RM 프로젝트 보안 표준에 따라 설치 과정을 로깅합니다.
```bash
# 설치 완료 로그 기록
echo "$(date): Bitwarden 서버 설치 완료 - 계정: bitwarden, 경로: /opt/ai4rm/docker/bitwarden" >> /var/log/ai4rm/security-audit.log
```

### 권한 검증
- 설치 후 권한이 올바르게 설정되었는지 확인합니다.
```bash
# 디렉토리 권한 확인
ls -la /opt/ai4rm/docker/bitwarden/
# CLI 실행 권한 확인  
ls -la /usr/local/bin/bw
```

## 문제해결

### 일반적인 문제
1. **권한 오류**: bitwarden 계정에 sudo 권한이 없는 경우
   - `/etc/sudoers.d/bitwarden-docker` 파일 설정 재확인
   
2. **포트 충돌**: 80/443 포트가 이미 사용 중인 경우
   - `netstat -tlnp | grep :80` 로 포트 사용 상태 확인
   
3. **메모리 부족**: 설치 중 메모리 부족 오류
   - `free -h` 로 메모리 상태 확인 후 스왑 공간 추가 고려

### 로그 확인
```bash
# Bitwarden 컨테이너 로그 확인
sudo docker compose -f /opt/ai4rm/docker/bitwarden/bwdata/docker/docker-compose.yml logs
```

## 참고자료
- [Bitwarden 공식 설치 가이드](https://bitwarden.com/help/install-on-premise-linux/)
- AI4RM 보안 정책: `/home/ben/projects/ai4rm/.github/copilot-instructions.md`

