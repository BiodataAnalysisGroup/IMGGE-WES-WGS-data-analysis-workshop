# Day 3

# Introduction to Common Workflow Language (CWL) and Docker

## Installation notes for WSL2 and Docker Desktop

To install the computational environment required for the execution of the present tutorial, follow the instructions of the [cwltool installation tutorial](https://github.com/common-workflow-language/cwltool), specifically the section on **MS Windows users**. There you will find detailed instructions on how to activate and install:

- Windows Subsystem for Linux 2 (WSL2)
- Docker Desktop
- Ubuntu or Debian from the Microsoft Store

## Overview

Reproducible research is an important part of good scientific practice. Establishing robust computational analysis workflow facilitates reproducible research. There has been a strong emphasis on establishing graphic workflow management systems, such as Galaxy and Knime. However, the landscape and complexity of Linux command line tools for workflow management is vast and complicated. We aim to take advantage of European/International leaders in the field to present current scientific workflow paradigms and a lead workshops to implement a basic workflow. Topics covered include workflows and containerization.

## Learning objectives

### Docker
- Understand what containerisation is, and why you might use it in bioinformatics.
- Be familiar with Docker; basic concepts and structure.
- Find and run containers built by other people.

### CWL
- Understand syntax and structure for CWL.
- Understand how to write CWL tool definitions for command line tools.
- Read and write CWL files written in YAML.
- Join CWL tools into a workflow.
- Use Docker with CWL to provide software dependencies and ensure reproducibility.

## Audience and requirements

This introductory tutorial is aimed towards bioinformaticians (graduate students and researchers), who are interested in becoming familiar to Docker based workflows in CWL.

## Prerequisites

- Experience in Shell; this includes basic commands (such as `ls`, `cp`, `mv`, `nano/vim`) and operations such as (`apt`, installing tools etc).

**Given that the workshop will offer hands-on training, participants will need to bring their own laptop with them.**
