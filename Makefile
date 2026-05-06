.PHONY: setup pipeline dashboard

VENV ?= .venv
ifeq ($(OS),Windows_NT)
PYTHON ?= $(VENV)/Scripts/python.exe
else
PYTHON ?= $(VENV)/bin/python
endif
PIP ?= $(PYTHON) -m pip

setup:
	python -m venv $(VENV)
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) analysis.py

dashboard:
	$(PYTHON) -m streamlit run dashboard.py --server.address 0.0.0.0 --server.port 8501
