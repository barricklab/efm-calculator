from django.conf.urls import patterns, url

from efm import views

urlpatterns = patterns('',
    url(r'^$', views.get_sequence, name='index'),
)
