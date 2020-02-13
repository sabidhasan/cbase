from django.shortcuts import render
from django.http import HttpResponse

def explore(request):
  return HttpResponse('This will list all chemicals.')

def my_chemicals(request):
  return HttpResponse('This will list MY chemicals!')
