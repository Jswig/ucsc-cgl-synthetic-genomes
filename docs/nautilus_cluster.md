# Nautilsu cluster usage

## Cluster

Getting data: using AWS S3 with Ceph file system
~15GB is reasonable in that
Write a small scipt to retrieve data from the S3 instance

## Kubernetes

Run as jobs, not pods - pods are limited to 2 cores and 8GB Ram

Only submit GPU requests under limits (not request).
Kill GPU jobs under 60% utilization.

- use `python -u` to unbuffer output
- use comet.ml for model monitoring

### Storage

- don't install AWS endpoints!
- just alias the URL for the endpoint
