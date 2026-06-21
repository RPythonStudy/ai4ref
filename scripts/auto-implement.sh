#!/bin/bash
# ai4ref auto-implement — tasks.md 의 미완료 묶음을 Issue-First 원칙대로 구현·이슈연결·종료 무인 반복.
# (cortex auto-implement.sh 포팅 — v0.11.3 하이픈 /speckit-implement)
#
# 사용법: ./scripts/auto-implement.sh specs/<NNN-feature>/tasks.md
# 예시:   ./scripts/auto-implement.sh specs/001-kq-pico-authoring/tasks.md
#
# ⚠️ 주의: --dangerously-skip-permissions(권한 우회·무인 자율) + 자동 commit/push + 이슈 생성/종료.
#         첫 실행은 반드시 지켜볼 것. 전제: gh 인증, GitHub 원격, CLAUDE.md Issue-First 원칙.

TASK_FILE="${1:?사용법: $0 <tasks.md 경로>}"
[ -f "$TASK_FILE" ] || { echo "❌ 파일 없음: $TASK_FILE"; exit 1; }

# 미완료 '- [ ]' 가 남아있는 동안 반복
while grep -Fq -e "- [ ]" "$TASK_FILE"; do
  TOTAL=$(grep -cE "^[[:space:]]*- \[[ xX]\]" "$TASK_FILE" || echo 0)
  REMAIN=$(grep -cF -e "- [ ]" "$TASK_FILE" || echo 0)
  DONE=$((TOTAL - REMAIN))
  PCT=0; [ "$TOTAL" -gt 0 ] && PCT=$((DONE * 100 / TOTAL))
  NEXT=$(grep -m1 -F -e "- [ ]" "$TASK_FILE" | sed -E 's/^[[:space:]]*- \[ \][[:space:]]*//')

  echo ""
  echo "======================================================="
  echo "📊 진행률: $DONE / $TOTAL ($PCT%)"
  echo "🎯 다음: $NEXT"
  echo "======================================================="
  echo "🚀 컨텍스트 초기화 후 다음 묶음 처리..."

  claude -p "/speckit-implement \"1. ${TASK_FILE} 에서 아직 완료되지 않은([ ]) 다음 미완료 묶음(같은 Phase 또는 User Story) 전체만 구현해줘(그 다음 묶음은 건드리지 마). 2. 테스트가 모두 통과하면 자동 커밋·푸시해. 3. CLAUDE.md 의 Issue-First 원칙대로 매칭 GitHub 이슈를 gh issue list 로 찾고, 없으면 gh issue create 로 생성한 뒤, 완료 시 gh issue close #번호 --reason completed 로 종료해. 4. 커밋 메시지 끝에 Closes #번호 를 붙이고 헌법 커밋 규칙(한글 서술·Co-Authored-By)을 따라. 5. **${TASK_FILE} 에서 완료한 항목을 반드시 [x] 로 체크해(루프 종료 조건).** 6. 완료된 Task ID·이슈 번호와 진행 상황을 짧게 보고해.\"" --dangerously-skip-permissions

  echo "⏳ 사이클 완료. 2초 후 남은 작업 재확인..."
  sleep 2
done

echo "🎉 모든 작업 완료!"
