FROM python:3.11-slim

# Prevent Python from writing .pyc files and enable unbuffered output
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PYTHONPATH="/app/src"

# Install system dependencies needed for RDKit, scientific libraries, and curl healthchecks
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    libxrender1 \
    libxext6 \
    libsm6 \
    libgl1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy dependency specifications first to leverage Docker layer caching
COPY requirements.txt pyproject.toml ./

# Install CPU PyTorch wheel directly with retry tolerance to avoid 1.5GB CUDA wheel download timeouts
RUN pip install --no-cache-dir --default-timeout=180 --retries 5 torch --index-url https://download.pytorch.org/whl/cpu

# Install remaining dependencies with timeout tolerance
RUN pip install --no-cache-dir --default-timeout=180 --retries 5 -r requirements.txt

# Copy application source code, configuration, artifacts, and data
COPY src/ /app/src/
COPY config/ /app/config/
COPY tests/ /app/tests/
COPY api/ /app/api/
COPY ui/ /app/ui/
COPY scripts/ /app/scripts/
COPY artifacts/ /app/artifacts/
COPY data/ /app/data/
COPY README.md DEMO_GUIDE.md /app/

# Create runtime directories (skills_cache is auto-seeded by VerifiedProceduralSkillLibrary on first startup)
RUN mkdir -p /app/skills_cache

# Install local package in editable mode
RUN pip install --no-cache-dir -e .

EXPOSE 8000 8501

CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
