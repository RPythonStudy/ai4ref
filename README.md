# 프로젝트 목적
AI를 활용하여 논문자동수집, 분석, 보고서 작성 등의 연구업무를 자동화하는 것을 목표로 합니다. 

# 프로젝트 개요
- 로컬참조는 github상태에서는 링크되지 않으므로 사용자분들은 원격참조를 사용하시기 바랍니다.
- AI Agent를 위한 로컬참조: [Project-Overview](wiki/Project-Overview.md)
- 사용자를 위한 원격참조: [Project Overview](https://github.com/RPythonStudy/rpy-quarto-template/wiki/Project-Overview)

# README.md 목적
- 이 프로젝트에서 `README.md`파일은 개발자가 AI agent에게 스크립트 자동생성을 요청할 때 AI agent가 참고해야 할 지침을 제공함을 목적으로 합니다.
- VS Code에서 copilot를 사용하는 경우를 위해서 프로젝트 폴더의 `.github/copilot-instructions.md`로부터 `README.md`로 심볼릭 링크를 구성하였습니다.

# Instructions

## Logging in R
- `src/R/logger.R` 사용
- .env파일로부터 PROJECT_NAME, LOG_LEVEL, SERVICE_LOG_PATH 로딩
- 함수: log_debug, log_info, log_warn, log_error
- log_fatal은 파이썬과 함수명을 통일하기 위해 log_critical로 래핑
- 정책상 logs 폴더 직접생성해야 하고, 권한 오류 발생가능
- 로컬참조: [R-Logging-Reference](wiki/R-Logging-Reference.md)
- 원격참조: [R Logging Reference](https://github.com/RPythonStudy/rpy-quarto-template/wiki/R-Logging-Reference)
- Home참조: [Wiki Home](https://github.com/RPythonStudy/rpy-quarto-template/wiki/Home)



## Logging in Python 
- 로깅을 구현한 파일 위치: src/common/logger.py
- syspath에 src를 추가했으므로 `from common.logger import log_info`와 같이 import
- .env파일로부터 PROJECT_NAME, LOG_LEVEL, SERVICE_LOG_PATH 로딩
- config/logging.yml 파일로부터 로그포맷과 같은 복잡한 정책 관리
- 개인정보호 logging framework를 차용해 왔기 때문에 정책상 logs 폴더 직접생성해야 함 (단 설치 자동화 스크립트에서는 자동 생성구현됨). 폴더의 위치에 따라서는 권한문제 발생 가능
- 구체적인 예시는 아래의 참조를 참고
- 참조: [Python-Logging-Reference](wiki/Python-Logging-Reference.md)


## 스크립트 제안 지침
- **극단적으로 간결·직관·디버그 친화적** 코드를 우선 제안
- Windows/Linux/Mac 간 경로 차이와 권한 이슈를 **사전 고지**하고, 필요 시 대체 경로 예시를 함께 제안
- 에러 메시지(예: 환경변수 치환 실패, 권한 오류)를 기반으로 **수정 가이드라인**까지 함께 제안

