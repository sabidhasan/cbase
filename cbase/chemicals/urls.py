from django.urls import path
from . import views

urlpatterns = [
  path('explore', views.explore, name='explore'),
  path('my-chemicals', views.my_chemicals, name='my-chemicals'),
]
