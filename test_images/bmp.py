from PIL import Image

def resize_and_convert_to_bmp(input_image_path, output_bmp_path, new_width, new_height):
    # 打開原始圖像
    image = Image.open(input_image_path)
    print(image.size)
    # 調整圖像大小
    resized_image = image.resize((new_width, new_height))
    
    # 將圖像轉換為 BMP 格式
    bmp_image = resized_image.convert("RGB")
    
    # 保存 BMP 圖像
    bmp_image.save(output_bmp_path)

    # 修改 BMP 圖像的 header
    with open(output_bmp_path, "r+b") as bmp_file:
        # 讀取文件頭部
        bmp_file.seek(18)  # 寬度信息在文件的 18-21 字節
        bmp_file.write(new_width.to_bytes(4, byteorder='little'))  # 將新的寬度寫入文件

        bmp_file.seek(22)  # 高度信息在文件的 22-25 字節
        bmp_file.write(new_height.to_bytes(4, byteorder='little'))  # 將新的高度寫入文件

def read_bmp_header(bmp_file_path):
    with open(bmp_file_path, "rb") as bmp_file:
        # 讀取 BMP 文件頭
        bmp_header = bmp_file.read(54)  # BMP 文件頭為固定 54 個字節

        # 提取寬度和高度信息
        width_bytes = bmp_header[18:22]  # 寬度信息在文件的 18-21 字節
        height_bytes = bmp_header[22:26]  # 高度信息在文件的 22-25 字節

        # 將字節數組轉換為整數
        width = int.from_bytes(width_bytes, byteorder='little')  # 寬度使用小端序
        height = int.from_bytes(height_bytes, byteorder='little')  # 高度使用小端序

        # 輸出尺寸相關訊息
        print("BMP 圖像尺寸信息:")
        print("寬度:", width, "像素")
        print("高度:", height, "像素")

input_image_path = "Anya.jpg"   # 輸入圖像路徑
output_bmp_path = "anya.bmp"  # 輸出 BMP 圖像路徑

# 這裡兩者相等是為了配合 bmp.cpp
new_width = 800   # 新的寬度
new_height = 800  # 新的高度

# 調整大小並轉換為 BMP 格式並修改 header
resize_and_convert_to_bmp(input_image_path, output_bmp_path, new_width, new_height)

# 讀取 BMP header 並輸出尺寸相關訊息
read_bmp_header(output_bmp_path)