@echo off
REM
REM Create temporary environment variables for JDK.
REM
REM This file is automatically generated. Please, do not edit this file.
REM


SETX JAVA_HOME %~dp0
for /f "skip=2 tokens=3*" %a in ('reg query HKCU\Environment /v PATH') do @if [%b]==[] ( @setx PATH "%~a;%JAVA_HOME\bin" ) else ( @setx PATH "%~a %~b;%JAVA_HOME\bin" )