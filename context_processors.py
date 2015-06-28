__author__ = 'ben'

from django.conf import settings

#######################################################################################
# This function is here specifically for tracking web traffic with google analytics.
# 'GOOGLE_ANALYTICS_PROPERTY_ID' and 'GOOGLE_ANALYTICS_DOMAIN' must be defined in
# django's 'settings.py' to enable functionality.
#######################################################################################

def google_analytics(request):
    """
    Use the variables returned in this function to
    render your Google Analytics tracking code template.
    """
    ga_prop_id = getattr(settings, 'GOOGLE_ANALYTICS_PROPERTY_ID', False)
    ga_domain = getattr(settings, 'GOOGLE_ANALYTICS_DOMAIN', False)
    print ga_prop_id
    if not settings.DEBUG and ga_prop_id and ga_domain:
        return {
            'GOOGLE_ANALYTICS_PROPERTY_ID': ga_prop_id,
            'GOOGLE_ANALYTICS_DOMAIN': ga_domain,
        }
    return {}
