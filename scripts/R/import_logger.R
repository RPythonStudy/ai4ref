# 리이브러리 일괄호출
library(here)

# 사용자정의함수 호출
invisible(lapply(
  list.files(here::here("src", "R"), pattern = "\\.R$", full.names = TRUE),
  source,
  encoding = "UTF-8"
))

# 로거 초기화
init_logger()

log_info("사용자정의함수 호출 완료")
