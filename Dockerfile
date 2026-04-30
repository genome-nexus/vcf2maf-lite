FROM python:3.10-slim

WORKDIR /vcf2maf_lite

RUN apt-get update && apt-get install -y curl \
    && curl -sSL https://install.python-poetry.org | python3 -
ENV PATH="/root/.local/bin:$PATH"
COPY pyproject.toml poetry.lock* ./
RUN poetry install

COPY . .
ENV PYTHONPATH=/vcf2maf_lite/vcf2maf_lite

ENTRYPOINT ["poetry", "run"]
CMD ["python3", "vcf2maf_lite/vcf2maf_lite.py"]


