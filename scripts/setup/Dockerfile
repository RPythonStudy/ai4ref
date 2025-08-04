FROM test-pg:latest

# 원하는 UID/GID로 변경 (예: 1022)
RUN groupmod -g 1022 postgres && \
    usermod -u 1022 -g 1022 postgres && \
    chown -R postgres:postgres /var/lib/postgresql

USER postgres
