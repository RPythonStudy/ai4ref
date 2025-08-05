# Zotero Web API Key 발급 및 `.env` 설정 가이드

Zotero Web API를 활용하여 논문 서지정보를 자동으로 저장하려면, **개인 API Key**가 필요합니다. 이 문서는 Zotero API Key 발급 과정과 `.env` 파일에 설정하는 방법을 안내합니다.

---

## 1. Zotero 웹사이트 접속 및 로그인

- 주소: https://www.zotero.org
- 계정이 없다면 무료 회원가입
- 로그인 후 다음 단계 진행

---

## 2. API 키 발급 페이지로 이동

- 상단 메뉴: `Settings` → `security` 탭  
  또는 바로 접속:  
  https://www.zotero.org/settings/security

---

## 3. 새 API 키 생성

- `Create new private key` 버튼 클릭

---

## 4. 권한 및 옵션 설정

- Name:  
  예: `ai4ref`

- Access:
  - `Allow library access` 체크

- Permissions:
  - `Read/Write` 선택

- Group Access:
  - 그룹 Zotero 사용 시 권한 부여 가능 (선택)

- 설정 완료 후 `Save Key` 클릭

---

## 5. API Key 복사 및 저장

- 생성된 키가 한 번만 표시되므로 반드시 복사해 둘 것 (Bitwarden 사용추천)

예시:
```
ZOTERO_API_KEY="dh1qR4KshQrUBPq7Jc33YW1s"
ZOTERO_USER_ID="17692286"
```

---

## 6. `.env` 파일에 저장

프로젝트 루트 디렉토리에 있는 `.env` 파일에 아래와 같이 저장합니다:

```env
# ZOTERO API 설정
ZOTERO_API_KEY="dh1qR4KshQrUBPq7Jc33YW1s"
ZOTERO_USER_ID="17692286"
```

> ⚠️ `.env` 파일은 Git 등 버전관리에서 **절대 공유되지 않도록 `.gitignore`에 포함**되어야 합니다.

---

## 7. Python 코드에서 사용하는 방법 예시

```python
import os
from dotenv import load_dotenv

load_dotenv()

ZOTERO_API_KEY = os.getenv("ZOTERO_API_KEY")
ZOTERO_USER_ID = os.getenv("ZOTERO_USER_ID")

headers = {
    "Authorization": f"Bearer {ZOTERO_API_KEY}",
    "Content-Type": "application/json"
}

api_url = f"https://api.zotero.org/users/{ZOTERO_USER_ID}/items"
```

---

## 8. 참고 자료

- Zotero Web API 공식 문서:  
  https://www.zotero.org/support/dev/web_api/v3/start

- API 요청 샘플 (논문 저장):
```http
POST https://api.zotero.org/users/1234567/items
Authorization: Bearer BvGf1bcA9kRCytXq7EaKoBsh
Content-Type: application/json
```

---

## 🔐 보안 주의

- `.env` 파일은 민감정보를 포함하므로 **절대 외부에 유출되지 않도록 주의**하세요.
- 팀 프로젝트에서는 Vault 등 비밀 관리 도구 사용을 고려하세요.

