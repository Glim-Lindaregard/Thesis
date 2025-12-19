# Run Training: 
```bash
PYTHON_PATH slider_isaac_lab/src/slider_isaac_lab/go_to_point/scripts/train.py \
--task Slider-GoToPoint-v0\
--headless\
--num_envs=512
```


# Run Validation
```bash
PYTHON_PATH slider_isaac_lab/src/slider_isaac_lab/go_to_point/scripts/agent_test.py \
--task Slider-GoToPoint-v0 \
--headless \
--num_envs=1 \
--enable_cameras \
--policy_path=<YOUR_POLICY_PATH>
```
