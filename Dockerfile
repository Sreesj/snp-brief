# Variant Atlas: docker build -t variant-atlas . && docker run -p 8501:8501 -e NEBIUS_API_KEY=... variant-atlas
FROM python:3.12-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
ENV SNP_BRIEF_USER_CACHE=/data/cache
VOLUME ["/data"]
EXPOSE 8501
CMD ["streamlit", "run", "streamlit_app.py", "--server.address=0.0.0.0", "--server.port=8501", "--browser.gatherUsageStats=false"]
