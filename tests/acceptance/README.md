# Acceptance Truth (검증 참값)

재작성(rewrite) 코드가 **반드시 재현해야 하는 동작의 참값**. 헌법 XIII(Validation Triad)·
Development Workflow "검증 참값 우선"의 근거 산출물 — 코드는 신규 재작성하되 이 참값을
통과해야 비로소 기존 MVP 와 등가로 인정한다.

## 파일

- **`recall_baseline.yml`** — 검색 recall 참값(두-검증 A/B) + 파생 검색식 스냅샷.
  `alert.validate_kq` 가 이 점수를 재현해야 등가. 측정 2026-06-22: **검증 A 8/9 · B 2/2**
  (누락 19602972 = ICU 간접근거, 의도된 천장).
- **`gate_baseline.yml`** — 게이트 precision 참값. 재작성 게이트가 `landmark_true` 는
  통과(is_landmark=true), `reject_false_positive` 는 거부(false)해야 등가. 백필 FP 0.

## 검증 (등가 판정)

```bash
uv run python -m alert.validate_kq     # recall_baseline 과 대조 → A 8/9 · B 2/2
# 게이트: 재작성 후 landmark_true / reject_false_positive 케이스로 확인 (claude CLI·Haiku)
```

> ⚠️ 이 파일들은 **참값**이라 임의 갱신 금지. PICO OR-리스트 정련으로 recall 이 오르면
> (예: 8/9→9/9) 그때 의도적으로 갱신하고 사유를 커밋 메시지에 남긴다. (정련 = 헌법 III)
