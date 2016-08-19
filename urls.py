from django.conf.urls import url

from efm import views

urlpatterns = [
    url(r'^', views.get_sequence, name='index'),
]
