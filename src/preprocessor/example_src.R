# example.R
# 교육용 예시 사용자정의 함수들

#' 간단한 인사 함수
#' @param name 이름
#' @return 인사말 문자열
say_hello <- function(name = "World") {
  log_info(paste("say_hello 함수 호출됨, name:", name))
  return(paste("안녕하세요,", name, "님!"))
}

#' 두 수의 합을 계산하는 함수
#' @param a 첫 번째 수
#' @param b 두 번째 수
#' @return 두 수의 합
add_numbers <- function(a, b) {
  result <- a + b
  log_debug(paste("add_numbers:", a, "+", b, "=", result))
  return(result)
}

#' 데이터프레임을 요약하는 함수
#' @param df 데이터프레임
#' @return 요약 정보
summarize_data <- function(df) {
  log_info("데이터프레임 요약 시작")
  
  summary_info <- list(
    rows = nrow(df),
    cols = ncol(df),
    col_names = names(df),
    na_count = sum(is.na(df))
  )
  
  log_info(paste("데이터 요약 완료 - 행:", summary_info$rows, "열:", summary_info$cols))
  return(summary_info)
}

# 함수 로딩 완료 로그
log_info("example.R 함수들이 로딩되었습니다")
